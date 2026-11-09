source('R/01_helpers.R')
cat('WP0-WP2: VCF verification, haploid recoding, sample/site QC\n')
A<-read_gt_vcf(FILE_ALL); B<-read_gt_vcf(FILE_2020)
assert(ncol(A$G)==EXPECTED_ALL && ncol(B$G)==EXPECTED_2020,'Unexpected sample count')
assert(nrow(A$G)==EXPECTED_PASS_SITES && nrow(B$G)==EXPECTED_PASS_SITES,'Unexpected site count')
assert(identical(A$sites, B$sites),'Full and 2020 cohort site metadata differ')
# Confirm derived VCF REF alleles against the locally downloaded 3D7 chromosome-5 FASTA.
fa <- readLines(required(FASTA), warn=FALSE)
h <- grep('^>',fa)
target <- which(grepl(paste0('^>',REGION_CHROM,'([[:space:]]|$)'),fa[h]))
assert(length(target)==1L,'Chromosome 5 not uniquely found in FASTA header')
fstart<-h[target]+1L; fend<-if(target<length(h))h[target+1L]-1L else length(fa)
reference_chr<-paste0(fa[fstart:fend],collapse='')
obs_ref<-toupper(mapply(function(p,r)substr(reference_chr,p,p+nchar(r)-1L),
                         A$sites$pos,A$sites$ref,USE.NAMES=FALSE))
assert(all(obs_ref==toupper(A$sites$ref)),'Derived VCF REF does not match the downloaded 3D7 reference')
rm(fa,reference_chr)
assert(all(B$ids %in% A$ids),'2020 sample list is not strict subset of 282')
assert(identical(unname(A$G[,B$ids,drop=FALSE]),unname(B$G)),'2020 genotype values differ from full 282 data')
stateA<-state_from_lists(A$ids,'all'); stateB<-state_from_lists(B$ids,'2020')
metaA<-load_sample_metadata(A$ids,stateA); metaB<-load_sample_metadata(B$ids,stateB)
assert(all(metaB$study=='1306-PF-NG-NGWA-SM' & metaB$year==2020), '2020 provenance check failed')
expectedA<-c(Oyo=71,Ogun=38,Lagos=129,Osun=15,Ondo=29)
expectedB<-c(Oyo=16,Ogun=38,Lagos=21,Osun=15,Ondo=29)
assert(identical(as.integer(table(factor(metaA$state,STATES))),as.integer(expectedA)), 'All cohort state counts differ')
assert(identical(as.integer(table(factor(metaB$state,STATES))),as.integer(expectedB)), '2020 state counts differ')
sha <- function(p) { s<-system2('sha256sum',p,stdout=TRUE,stderr=FALSE); strsplit(trimws(s),' +')[[1]][1] }
manifest<-data.frame(role=c('282_PASS_VCF','119_PASS_VCF','sample_metadata','Fws'),
  path=c(FILE_ALL,FILE_2020,META,FWS), stringsAsFactors=FALSE)
manifest$sha256<-vapply(manifest$path,sha,character(1))
write_tsv(manifest,'input_manifest_SHA256.tsv','audit')
config_audit<-file.path(OUT,'audit','ANALYSIS_CONFIG_at_first_run.R')
if(!file.exists(config_audit))file.copy('R/00_config.R',config_audit,overwrite=FALSE)
log_mat <- function(X,met,label) {
  sample_missing<-colMeans(is.na(X$G))
  sample_ambig<-colMeans(X$ambiguous)
  sample_partial<-colMeans(X$partial)
  met$pass_missingness_haploid<-sample_missing
  met$pass_ambiguity_fraction<-sample_ambig
  met$pass_partial_fraction<-sample_partial
  met$callable_sites<-colSums(!is.na(X$G))
  met$cohort<-label
  met
}
metaA<-log_mat(A,metaA,'all282'); metaB<-log_mat(B,metaB,'primary2020')
siteA<-site_stats(A$G,metaA$state,'all282'); siteB<-site_stats(B$G,metaB$state,'primary2020')
siteA$pos<-A$sites$pos;siteB$pos<-B$sites$pos
siteA$ref<-A$sites$ref;siteB$ref<-B$sites$ref
siteA$alt<-A$sites$alt;siteB$alt<-B$sites$alt
write_tsv(metaA,'Table_S1_full_sample_QC.tsv')
write_tsv(metaB,'Table_S2_2020_sample_QC.tsv')
write_tsv(siteA,'Table_S3_full_site_QC.tsv')
write_tsv(siteB,'Table_S4_2020_site_QC.tsv')
cohorts<-list(all282=list(X=A,meta=metaA,site=siteA),primary2020=list(X=B,meta=metaB,site=siteB))
saveRDS(cohorts,file.path(OUT,'rds','QC_cohorts.rds'),compress='xz')
make_state_table<-function(m,cohort) {
  z<-do.call(rbind,lapply(STATES,function(s) {
    x<-m[m$state==s,,drop=FALSE]
    data.frame(cohort=cohort,state=s,n=nrow(x),
      fws_median=median(x$Fws),fws_q1=unname(quantile(x$Fws,.25)),fws_q3=unname(quantile(x$Fws,.75)),
      fws_ge095=sum(x$Fws>=.95),
      pass_missingness_median=median(x$pass_missingness_haploid),
      missingness_le10=sum(x$pass_missingness_haploid<=.10),
      missingness_gt20=sum(x$pass_missingness_haploid>.20),
      ambiguous_two_slot_calls=sum(x$pass_ambiguity_fraction*EXPECTED_PASS_SITES),
      studies=length(unique(x$study)),years=length(unique(na.omit(x$year))))
  }))
  z
}
t1<-rbind(make_state_table(metaA,'all282'),make_state_table(metaB,'primary2020'))
write_tsv(t1,'Table_1_study_population_and_QC.tsv')
filter_counts<-function(z,cohort) {
  ss<-z$site
  data.frame(cohort=cohort,phase=c('Pf8 global PASS','Observed >=2 allele codes','Biallelic','Biallelic MAC >=3',
            'Call rate >=90%','Per-state call >=80%','All prespecified criteria'),
    sites=c(nrow(ss),sum(ss$allele_count>=2),sum(ss$allele_count==2),
       sum(ss$allele_count==2 & !is.na(ss$minor_count) & ss$minor_count>=3),
       sum(ss$callable>=.90),
       sum(rowSums(ss[,STATES]>=.80)==length(STATES)),sum(ss$informative)))
}
write_tsv(rbind(filter_counts(cohorts$all282,'all282'),filter_counts(cohorts$primary2020,'primary2020')),
 'Table_2_site_attrition.tsv')
cat('FULL 282:  sample and site QC complete. Polymorphic after haploid consensus:',sum(siteA$allele_count>=2),'\n')
cat('PRIMARY 119: sample and site QC complete. Polymorphic after haploid consensus:',sum(siteB$allele_count>=2),'\n')
cat('PRIMARY 119: prespecified high-quality biallelic sites:',sum(siteB$informative),'\n')
if(requireNamespace('ggplot2',quietly=TRUE)) {
  library(ggplot2)
  sampleplot<-ggplot(t1,aes(x=state,y=n,fill=cohort))+geom_col(position='dodge',width=.72)+
    geom_text(aes(label=n),position=position_dodge(width=.72),vjust=-.25,size=3)+
    scale_fill_manual(values=c(all282='#3C5672',primary2020='#C07146'))+
    labs(title='PI4K study populations across Southwestern Nigeria',x=NULL,y='Samples',fill='Analysis population')+
    clean_theme()+scale_y_continuous(expand=expansion(mult=c(0,.15)))
  export_plot(sampleplot,'Figure_1_sampling_cohorts',8,4.8)
  qcp<-ggplot(metaB,aes(x=state,y=pass_missingness_haploid,color=state))+
    geom_boxplot(outlier.shape=NA,width=.50,color='#3B4551',fill='white')+
    geom_jitter(width=.11,height=0,alpha=.72,size=1.6)+
    geom_hline(yintercept=.10,linetype='dashed',color='#AA574D')+
    scale_y_continuous(labels=function(x) paste0(round(100*x),'%'))+
    labs(title='PI4K call missingness in the controlled 2020 cohort',
      subtitle='Dashed line: exploratory 10% missingness criterion',x=NULL,y='PASS-site missingness')+
    guides(color='none')+clean_theme()
  export_plot(qcp,'Figure_2_2020_sample_missingness',8,5)
  fp<-ggplot(metaB,aes(x=state,y=Fws,color=state))+geom_boxplot(outlier.shape=NA,
      color='#3B4551',fill='white',width=.5)+geom_jitter(width=.12,size=1.6,alpha=.75)+
      geom_hline(yintercept=.95,linetype='dashed',color='#AA574D')+
      guides(color='none')+labs(title='Within-host diversity metric (Fws), 2020 samples',
      subtitle='Fws >=0.95 is an exploratory low-diversity stratum, not proof of monoclonality',x=NULL,y='Fws')+clean_theme()
  export_plot(fp,'Figure_3_2020_Fws',8,5)
  lp<-ggplot(siteB,aes(x=pos,y=callable,color=informative))+
    geom_point(alpha=.5,size=.75)+geom_hline(yintercept=.90,linetype='dotted')+
    scale_color_manual(values=c('FALSE'='#B3B5BA','TRUE'='#AF5A3F'))+
    labs(title='PI4K PASS SNP callability across the 2020 cohort',
      subtitle='Highlighted sites meet prespecified descriptive criteria; not yet a barcode',
      x='3D7 chromosome 5 position (bp)',y='Fraction callable',color='Informative')+clean_theme()
  export_plot(lp,'Figure_4_PI4K_site_callability',9,4.7)
}
writeLines(capture.output(sessionInfo()),file.path(OUT,'audit','R_sessionInfo.txt'))
cat('WP0-WP2 COMPLETE — outputs:',OUT,'\n')
