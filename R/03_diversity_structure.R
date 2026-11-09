source('R/01_helpers.R')
cat('WP3: Allele frequencies, sample-size-aware diversity, exploratory PCoA and LD\n')
cohorts<-readRDS(required(file.path(OUT,'rds','QC_cohorts.rds')))
set.seed(SETTINGS$seed)
all_frequency<-list();diversity<-list(); rarefy<-list()
for (nm in names(cohorts)) {
  z<-cohorts[[nm]]; G<-z$X$G; m<-z$meta; ss<-z$site
  af<-list();ds<-list();rf<-list()
  for (s in STATES) {
    sub<-G[,m$state==s,drop=FALSE]
    # Global, NOT state-local, allele spectrum establishes whether ALT1 (code 1)
    # is the sole nonreference allele. Do not drop sites that are fixed (0 or 1)
    # within a particular state: these may be geographically informative.
    for (i in which(ss$allele_count==2L)) {
      global_codes <- sort(unique(as.integer(na.omit(G[i,]))))
      if (!identical(global_codes, c(0L,1L))) next
      x <- sub[i, ]
      n <- sum(!is.na(x))
      k <- sum(x==1L, na.rm=TRUE)
      ci <- if (n > 0L) wilson_ci(k,n) else c(NA_real_,NA_real_)
      af[[length(af)+1L]] <- data.frame(cohort=nm,state=s,pos=ss$pos[i],
        n_callable=n,n_alt=k,alt_frequency=if(n>0L) k/n else NA_real_,
        ci_low=ci[1],ci_high=ci[2],
        minor_allele_count=if(n>0L) min(k,n-k) else NA_integer_)
    }
    idx<-which(ss$informative)
    H<-vapply(idx,function(i) {
       x<-na.omit(sub[i,]); if (length(x)<2L) return(NA_real_)
       f<-as.numeric(table(x))/length(x); 1-sum(f*f)
    },numeric(1))
    ds[[length(ds)+1]]<-data.frame(cohort=nm,state=s,n_samples=ncol(sub),
      informative_sites=length(idx),number_callable_site_estimates=sum(is.finite(H)),
      mean_unbiased_site_heterozygosity=NA_real_) # estimated from observed haploid calls below
    # Unbiased haploid within-state allelic diversity n/(n-1)*(1-sum(p^2)), per site.
    div<-vapply(idx,function(i){x<-na.omit(sub[i,]);n<-length(x)
      if(n<2)return(NA_real_)
      p<-as.numeric(table(x))/n
      n/(n-1)*(1-sum(p^2))
    },numeric(1))
    ds[[length(ds)]]$mean_unbiased_site_heterozygosity<-mean(div,na.rm=TRUE)
    ds[[length(ds)]]$sum_site_diversity<-sum(div,na.rm=TRUE)
    ds[[length(ds)]]$n_nonmissing_diversity_sites<-sum(is.finite(div))
    ds[[length(ds)]]$within_state_variable_sites<-sum(vapply(idx, function(i) length(unique(na.omit(sub[i,])))>=2L,logical(1)))
    # Rarefaction, sampled without replacement, fixed n=min state sample size.
    target<-if(nm=='primary2020')15L else 15L
    vals<-rep(NA_integer_,200)
    if (ncol(sub)>=target) for (rep in seq_along(vals)) {
      smp<-sample(seq_len(ncol(sub)),target,replace=FALSE)
      vals[rep]<-sum(vapply(idx,function(i) length(unique(na.omit(sub[i,smp])))>=2L,logical(1)))
    }
    rf[[length(rf)+1]]<-data.frame(cohort=nm,state=s,n_rarefied=target,
      expected_variable_sites=mean(vals,na.rm=TRUE),
      lower_bootstrap=unname(quantile(vals,.025,na.rm=TRUE)),
      upper_bootstrap=unname(quantile(vals,.975,na.rm=TRUE)))
  }
  all_frequency[[nm]]<-if(length(af))do.call(rbind,af) else data.frame()
  diversity[[nm]]<-do.call(rbind,ds)
  rarefy[[nm]]<-do.call(rbind,rf)
}
AF<-do.call(rbind,all_frequency);DIV<-do.call(rbind,diversity);RARE<-do.call(rbind,rarefy)
write_tsv(AF,'Table_S5_per_state_allele_frequencies_with_Wilson_CI.tsv')
write_tsv(DIV,'Table_3_within_PI4K_allelic_diversity.tsv')
write_tsv(RARE,'Table_S6_rarefied_variable_site_richness.tsv')
# PCoA using unambiguous haploid calls and per-pair overlap; NOT whole-genome population structure.
Z<-cohorts$primary2020; ix<-which(Z$site$informative)
X<-Z$X$G[ix,,drop=FALSE]; ids<-colnames(X)
n<-ncol(X); dis<-matrix(0,n,n,dimnames=list(ids,ids));overlap<-matrix(0L,n,n)
for(i in seq_len(n-1L)) for(j in (i+1L):n) {
  keep<-!is.na(X[,i]) & !is.na(X[,j]);nn<-sum(keep)
  overlap[i,j]<-overlap[j,i]<-nn
  dis[i,j]<-dis[j,i]<-if(nn>=max(10L,ceiling(.2*length(ix))))mean(X[keep,i]!=X[keep,j]) else NA_real_
}
# Inadequate overlap means no robust ordination, rather than silently replacing missing distances.
if (length(ix)>=10L && !anyNA(dis)) {
  ord<-cmdscale(as.dist(dis),k=2,eig=TRUE,add=TRUE)
  P<-data.frame(sample_id=ids,state=Z$meta$state,axis1=ord$points[,1],axis2=ord$points[,2],
    Fws=Z$meta$Fws,missingness=Z$meta$pass_missingness_haploid)
  write_tsv(P,'Table_S7_PI4K_PCoA_coordinates.tsv')
  write_tsv(data.frame(sample_a=rep(ids,each=n),sample_b=rep(ids,n),distance=as.vector(dis),
    shared_sites=as.vector(overlap)), 'Table_S8_PCoA_pairwise_distance.tsv')
  if(requireNamespace('ggplot2',quietly=TRUE)) {
    library(ggplot2)
    p<-ggplot(P,aes(axis1,axis2,color=state))+geom_point(size=2.3,alpha=.85)+
      labs(title='PI4K-local PCoA in common-study 2020 cohort',
           subtitle='Unambiguous consensus calls; axes are NOT genome-wide ancestry',
           x='PCoA axis 1',y='PCoA axis 2')+clean_theme()
    export_plot(p,'Figure_5_2020_PI4K_PCoA',7,6)
  }
} else cat('PCoA NOT PRODUCED: insufficient informative sites or pairwise callable overlap\n')
# r2 within one linked locus, not independent-marker LD.
score<-Z$site$minor_count; score[!Z$site$informative]<-NA
ix<-head(order(-score, Z$site$pos, na.last=NA),20L)
if (length(ix)>=2L) {
  ld<-matrix(NA_real_,length(ix),length(ix),dimnames=list(as.character(Z$site$pos[ix]),as.character(Z$site$pos[ix])))
  for (i in seq_along(ix)) for(j in seq_along(ix)) {
    a<-Z$X$G[ix[i],];b<-Z$X$G[ix[j],];ok<-!is.na(a)&!is.na(b)
    if(sum(ok)>=20L && length(unique(a[ok]))==2L && length(unique(b[ok]))==2L)
      ld[i,j]<-suppressWarnings(cor(a[ok],b[ok])^2)
  }
  write_tsv(data.frame(pos=rownames(ld),ld,check.names=FALSE),'Table_S9_PI4K_LD_r2.tsv')
  if(requireNamespace('ggplot2',quietly=TRUE)) {
    df<-as.data.frame(as.table(ld));names(df)<-c('site_a','site_b','r2')
    library(ggplot2)
    p<-ggplot(df,aes(site_a,site_b,fill=r2))+geom_tile()+coord_equal()+
      scale_fill_gradient(low='white',high='#274D72',na.value='#DDDDDD',limits=c(0,1))+
      labs(title='PI4K linkage disequilibrium (r²)',subtitle='Exploratory top-20 common variants; single linked gene',
        x='Position (bp)',y='Position (bp)',fill='r²')+
      theme_minimal(base_size=10)+theme(axis.text.x=element_text(angle=55,hjust=1),panel.grid=element_blank())
    export_plot(p,'Figure_6_PI4K_LD',8,7)
  }
}
# Exploratory complete-case mini-haplotypes, never phase inference from mixed infections.
ix<-head(order(-Z$site$minor_count,na.last=NA),8L)
ix<-ix[Z$site$informative[ix]]
if(length(ix)>=3L) {
  H<-Z$X$G[ix,,drop=FALSE]
  cc<-colSums(is.na(H))==0L
  if(sum(cc)>=20L) {
    h<-apply(H[,cc,drop=FALSE],2,paste,collapse=':')
    hap<-data.frame(sample_id=colnames(H)[cc],state=Z$meta$state[cc],haplotype_pattern=h)
    write_tsv(hap,'Table_S10_exploratory_PI4K_minihaplotypes.tsv')
  }
}
cat('WP3 COMPLETE — allele frequency, diversity, LD and conditional ordination outputs in', OUT,'\n')
