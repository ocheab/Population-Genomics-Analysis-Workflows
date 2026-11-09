source('R/01_helpers.R')
cat('WP4: Hudson-type haploid FST, isolate-level bootstrap, label permutation and site associations\n')
cohorts<-readRDS(required(file.path(OUT,'rds','QC_cohorts.rds')))
set.seed(SETTINGS$seed+4L)
hudson <- function(A,B) {
  assert(nrow(A)==nrow(B),'FST position mismatch')
  nA<-rowSums(!is.na(A));nB<-rowSums(!is.na(B))
  pA<-rowSums(A==1L,na.rm=TRUE)/pmax(1,nA)
  pB<-rowSums(B==1L,na.rm=TRUE)/pmax(1,nB)
  codesok<-rowSums(!is.na(A) & A!=0L & A!=1L)==0L &
           rowSums(!is.na(B) & B!=0L & B!=1L)==0L
  denom<-pA*(1-pB)+pB*(1-pA)
  use<-nA>=2L & nB>=2L & denom>0 & codesok
  if(sum(use)<2L) return(c(fst=NA_real_,sites=sum(use)))
  num<-(pA[use]-pB[use])^2-pA[use]*(1-pA[use])/(nA[use]-1)-pB[use]*(1-pB[use])/(nB[use]-1)
  c(fst=sum(num)/sum(denom[use]),sites=sum(use))
}
pair_fst<-function(G,state,cohort,bootstrap=0L,permutations=0L) {
  out<-list();pairs<-combn(STATES,2,simplify=FALSE)
  for (pair in pairs) {
    ids1<-which(state==pair[1]); ids2<-which(state==pair[2]);
    if(length(ids1)<4 || length(ids2)<4)next
    U<-G[,ids1,drop=FALSE];V<-G[,ids2,drop=FALSE]
    val<-hudson(U,V)
    bt<-numeric(0); pm<-numeric(0)
    if(bootstrap>0) bt<-replicate(bootstrap, {
      a<-U[,sample(seq_len(ncol(U)),ncol(U),replace=TRUE),drop=FALSE]
      b<-V[,sample(seq_len(ncol(V)),ncol(V),replace=TRUE),drop=FALSE]
      hudson(a,b)['fst']
    })
    if(permutations>0) {
      W<-cbind(U,V);n1<-ncol(U)
      pm<-replicate(permutations, {
        ix<-sample(seq_len(ncol(W)))
        hudson(W[,ix[seq_len(n1)],drop=FALSE],W[,ix[-seq_len(n1)],drop=FALSE])['fst']
      })
    }
    out[[length(out)+1]]<-data.frame(cohort=cohort,state_a=pair[1],state_b=pair[2],
      n_a=ncol(U),n_b=ncol(V),n_informative_pair_sites=unname(val['sites']),
      hudson_fst=unname(val['fst']),
      bootstrap_ci_low=if(length(bt)&&sum(is.finite(bt))>20)unname(quantile(bt,.025,na.rm=TRUE)) else NA_real_,
      bootstrap_ci_high=if(length(bt)&&sum(is.finite(bt))>20)unname(quantile(bt,.975,na.rm=TRUE)) else NA_real_,
      two_sided_permutation_p=if(length(pm)&&is.finite(val['fst']))(1+sum(abs(pm)>=abs(val['fst']),na.rm=TRUE))/(1+sum(is.finite(pm))) else NA_real_)
  }
  do.call(rbind,out)
}
Z<-cohorts$primary2020
# Use ALL observed biallelic PASS positions without ranking geographic SNP association first.
idx<-which(Z$site$allele_count==2L)
G<-Z$X$G[idx,,drop=FALSE]
FST<-pair_fst(G,Z$meta$state,'primary2020',SETTINGS$bootstrap_reps,SETTINGS$fst_perm_reps)
FST$BH_q<-p.adjust(FST$two_sided_permutation_p,method='BH')
write_tsv(FST,'Table_4_primary2020_pairwise_Fst.tsv')
# sensitivity: no pairwise permutation testing due small state counts, effect sizes only.
scenario<-list(
 missingness10=Z$meta$pass_missingness_haploid<=.10,
 fws95=Z$meta$Fws>=.95,
 joint=Z$meta$pass_missingness_haploid<=.10 & Z$meta$Fws>=.95)
sens<-list()
for(nm in names(scenario)) {
  ix_s<-which(scenario[[nm]])
  a<-pair_fst(G[,ix_s,drop=FALSE],Z$meta$state[ix_s],nm,bootstrap=150L,permutations=0L)
  if(!is.null(a))sens[[nm]]<-a
}
A<-cohorts$all282
idx_all<-which(A$site$allele_count==2L)
sens$all282<-pair_fst(A$X$G[idx_all,,drop=FALSE],A$meta$state,'all282_uncontrolled',bootstrap=150L)
if(length(sens))write_tsv(do.call(rbind,sens),'Table_S11_Fst_sensitivity.tsv')
# Whole-locus global genotype-centroid R2, within the common 2020 study-year.
# Scores are descriptive; missing values imputed by the POOLED variant mean only for this test.
R2_global <- function(X,state) {
  mu<-colMeans(X);TSS<-sum(sweep(X,2,mu,'-')^2)
  if(TSS==0)return(NA_real_)
  BTS<-0
  for(s in STATES) {xx<-X[state==s,,drop=FALSE];if(nrow(xx))BTS<-BTS+nrow(xx)*sum((colMeans(xx)-mu)^2)}
  BTS/TSS
}
if(nrow(G)>0L) {
  X<-t(G)
  # binary 0/1 biallelic only; pooled-mean imputation prevents forced ALT states.
  for(j in seq_len(ncol(X))) {miss<-is.na(X[,j]);if(any(miss)) {
      obs<-X[!miss,j]; X[miss,j]<-if(length(obs))mean(obs) else 0
  }}
  obsR2<-R2_global(X,Z$meta$state)
  nullR2<-replicate(SETTINGS$global_perm_reps,R2_global(X,sample(Z$meta$state)))
  global<-data.frame(cohort='primary2020',method='genotype_centroid_R2_pooled_missing_imputation',
    n_samples=nrow(X),n_snps=ncol(X),R2=obsR2,
    permutations=SETTINGS$global_perm_reps,
    permutation_p=(1+sum(nullR2>=obsR2,na.rm=TRUE))/(1+sum(is.finite(nullR2))))
  write_tsv(global,'Table_S12_global_PI4K_variance_permutation.tsv')
}
# Site-level association: 5 x 2 per-ISOLATE allele table, exact Monte Carlo, BH across all sites tested.
assoc<-list()
for(i in which(Z$site$allele_count==2L)) {
  x<-Z$X$G[i,];u<-sort(unique(na.omit(x)))
  if(!identical(as.integer(u),c(0L,1L)))next
  valid<-!is.na(x);tab<-table(factor(Z$meta$state[valid],levels=STATES),factor(x[valid],levels=c(0L,1L)))
  if(any(rowSums(tab)<5L) || min(colSums(tab))<3L)next
  # 2000 Monte Carlo replicates, not an asymptotic small-cell test.
  test<-tryCatch(suppressWarnings(fisher.test(tab,simulate.p.value=TRUE,B=1999)),error=function(e)NULL)
  if(is.null(test))next
  assoc[[length(assoc)+1]]<-data.frame(pos=Z$site$pos[i],n_callable=sum(valid),
    n_alt=sum(x==1L,na.rm=TRUE),p_montecarlo=unname(test$p.value))
}
if(length(assoc)) {
  AS<-do.call(rbind,assoc);AS$BH_q<-p.adjust(AS$p_montecarlo,'BH')
  AS<-AS[order(AS$BH_q,AS$p_montecarlo,AS$pos),]
  write_tsv(AS,'Table_S13_exploratory_site_association_BH.tsv')
}
if(requireNamespace('ggplot2',quietly=TRUE)) {
  library(ggplot2)
  q<-FST
  q$contrast<-factor(paste(q$state_a,q$state_b,sep=' vs '),
      levels=rev(paste(q$state_a,q$state_b,sep=' vs ')))
  p<-ggplot(q,aes(hudson_fst,contrast))+geom_vline(xintercept=0,linetype='dashed',color='grey45')+
    geom_errorbarh(aes(xmin=bootstrap_ci_low,xmax=bootstrap_ci_high),height=.20,color='#52708B')+
    geom_point(size=2.5,color='#BC5F45')+
    labs(title='Pairwise PI4K Hudson-type FST, 2020 common-study samples',
      subtitle='95% isolate-bootstrap intervals; linked SNPs are not independent units',
      x='FST (signed estimates retained)',y=NULL)+clean_theme()
  export_plot(p,'Figure_7_pairwise_PI4K_Fst',9,6)
}
cat('WP4 COMPLETE — FST, permutation and site-association outputs in',OUT,'\n')
