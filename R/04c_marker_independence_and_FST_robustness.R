# STEP 26 -- local, reproducible marker-pattern, LD and differentiation robustness review.
# Run from /mnt/c/PI4K_Nigeria_Study/PI4K_Southwest_R_Analysis_Suite
# Uses only audited RDS and Step25 tables. No source file alteration or download.
source('R/01_helpers.R')
cat('STEP 26: four-marker dependence, pattern completeness, FST sensitivity audit (OFFLINE)\n')
Z <- readRDS(required(file.path(OUT, 'rds', 'QC_cohorts.rds')))$primary2020
Y <- as.character(Z$meta$state)
ID <- as.character(Z$meta$sample_id)
assert(length(Y)==119L && length(ID)==119L, 'Expected 119 primary2020 isolates')
assert(!anyDuplicated(ID), 'Duplicate sample IDs')
assert(all(Y %in% STATES), 'Unknown state in 2020 cohort')
MARKERS <- c(409688L, 410338L, 410817L, 413463L)
ix <- match(MARKERS, Z$site$pos)
assert(!anyNA(ix), 'Missing required marker coordinate')
assert(all(Z$site$informative[ix]), 'One or more four-marker sites no longer meet frozen QC policy')
G <- t(Z$X$G[ix, , drop=FALSE]) # isolates x marker positions
storage.mode(G) <- 'integer'
colnames(G) <- as.character(MARKERS)
rownames(G) <- ID
assert(all(is.na(G) | G %in% c(0L, 1L)), '4-marker genotypes require additional allele recoding')
assert(identical(colnames(G),as.character(MARKERS)), 'Unexpected marker ordering')

# Data completeness at four selected sites; not a substitute for PI4K-wide missingness.
n_called <- rowSums(!is.na(G))
complete <- n_called == ncol(G)
pattern <- apply(G, 1, function(z) if(anyNA(z)) NA_character_ else paste(z, collapse='-'))
pattern_df <- data.frame(sample_id=ID, state=Y, n_markers_called=n_called,
   four_marker_complete=complete, four_marker_consensus_pattern=pattern,
   Fws=Z$meta$Fws, PI4K_PASS_missingness=Z$meta$pass_missingness_haploid,
   stringsAsFactors=FALSE)
write_tsv(pattern_df, 'Table_S30_2020_four_marker_sample_patterns.tsv')

# Include explicit zero-frequency pattern/state combinations to avoid omission.
all_patterns <- sort(unique(pattern[!is.na(pattern)]))
pat_state <- expand.grid(state=STATES, pattern=all_patterns, stringsAsFactors=FALSE)
pat_state$n <- mapply(function(st, pa) sum(Y==st & pattern==pa,na.rm=TRUE), pat_state$state, pat_state$pattern)
pat_state$state <- factor(pat_state$state, levels=STATES)
state_complete_n <- vapply(STATES,function(s)sum(Y==s & complete),integer(1))
pat_state$n_state_complete <- state_complete_n[as.character(pat_state$state)]
pat_state$proportion_among_complete <- ifelse(pat_state$n_state_complete>0,pat_state$n/pat_state$n_state_complete,NA_real_)
write_tsv(pat_state, 'Table_S31_2020_four_marker_consensus_patterns_by_state.tsv')

call_df <- do.call(rbind,lapply(STATES, function(s) {
  m <- Y==s; n<-sum(m); nc<-sum(complete[m]); ci<-wilson_ci(nc,n)
  data.frame(state=s, n_samples=n, n_complete_four_marker=nc,
   complete_call_fraction=nc/n, wilson_95_low=ci[1], wilson_95_high=ci[2],
   mean_n_markers_called=mean(n_called[m]),
   median_PI4K_PASS_missingness=median(Z$meta$pass_missingness_haploid[m]),
   median_Fws=median(Z$meta$Fws[m]), stringsAsFactors=FALSE)
}))
write_tsv(call_df, 'Table_S33_2020_four_marker_completeness_by_state.tsv')

# Per-pair linkage, including within-Ogun dependence of the two Ogun-enriched variants.
ld_rows <- list()
for (subset_label in c('All_2020', 'Ogun_2020')) {
  idx <- if(subset_label=='All_2020') seq_len(nrow(G)) else which(Y=='Ogun')
  for(p in combn(seq_along(MARKERS),2L,simplify=FALSE)) {
    a<-G[idx,p[1]]; b<-G[idx,p[2]]; ok<-!is.na(a)&!is.na(b)
    n<-sum(ok); variable_a<-length(unique(a[ok]))==2L; variable_b<-length(unique(b[ok]))==2L
    r2 <- if(n>=8L && variable_a && variable_b) suppressWarnings(cor(a[ok],b[ok])^2) else NA_real_
    # Co-occurrence counts for all allele combinations, regardless of variability.
    n00<-sum(a[ok]==0L & b[ok]==0L)
    n01<-sum(a[ok]==0L & b[ok]==1L)
    n10<-sum(a[ok]==1L & b[ok]==0L)
    n11<-sum(a[ok]==1L & b[ok]==1L)
    fisher_p<-if(variable_a && variable_b && n>0L)
      fisher.test(matrix(c(n00,n10,n01,n11),nrow=2L))$p.value else NA_real_
    ld_rows[[length(ld_rows)+1L]]<-data.frame(subset=subset_label,
      marker_a=MARKERS[p[1]],marker_b=MARKERS[p[2]],n_joint_called=n,
      n00=n00,n01=n01,n10=n10,n11=n11,
      r_squared=r2, fisher_exact_p_descriptive=fisher_p,
      note='Four markers in one linked PI4K region; no LD-independent replication',
      stringsAsFactors=FALSE)
  }
}
ld <- do.call(rbind,ld_rows)
write_tsv(ld,'Table_S32_2020_four_marker_pairwise_LD_and_cooccurrence.tsv')

a<-G[Y=='Ogun','409688']; b<-G[Y=='Ogun','413463'];ok<-!is.na(a)&!is.na(b)
cooc <- data.frame(first_marker_409688=c('A (code 0)','A (code 0)','T (code 1)','T (code 1)'),
 second_marker_413463=c('G (code 0)','A (code 1)','G (code 0)','A (code 1)'),
 n_Ogun=c(sum(a[ok]==0 & b[ok]==0),sum(a[ok]==0 & b[ok]==1),
          sum(a[ok]==1 & b[ok]==0),sum(a[ok]==1 & b[ok]==1)),
 n_joint_called=rep(sum(ok),4), stringsAsFactors=FALSE)
write_tsv(cooc,'Table_S35_Ogun_2020_two_marker_cooccurrence.tsv')

# Reconcile strict, broad-consensus and 282-site FST as distinct analyses.
ft <- read_tsv(file.path(OUT,'tables','Table_S26_PI4K_Fst_primary_and_sensitivity_AUDITED.tsv'))
required_columns<-c('cohort','marker_policy','state_a','state_b','hudson_Fst_signed',
 'n_pair_informative_sites','BH_q_10_pairs')
assert(all(required_columns %in% names(ft)), 'Step25 FST table has unexpected columns')
ft$scenario <- ifelse(ft$cohort=='2020_same_study' & ft$marker_policy=='strict',
 '2020_four_strict',ifelse(ft$cohort=='2020_same_study',
 '2020_consensus_43_exploratory','282_strict_18_confounded'))
key<-paste(ft$state_a,ft$state_b,sep='__')
expected_pairs <- apply(combn(STATES, 2L), 2L, paste,collapse='__')
assert(all(expected_pairs %in% key[ft$scenario=='2020_four_strict']),
 'FST primary result missing pair contrasts')
robust <- data.frame(state_a=character(),state_b=character(),
 scenario=character(),n_selected_sites=integer(),n_pair_informative_sites=integer(),
 Hudson_Fst=numeric(),BH_q=numeric(),CI_low=numeric(),CI_high=numeric(),
 interval_available=logical(),stringsAsFactors=FALSE)
for(i in seq_len(nrow(ft))) {
  low<-as.numeric(ft$isolate_bootstrap_CI_low[i]);hi<-as.numeric(ft$isolate_bootstrap_CI_high[i])
  robust<-rbind(robust,data.frame(state_a=ft$state_a[i],state_b=ft$state_b[i],
   scenario=ft$scenario[i],n_selected_sites=as.integer(ft$n_selected_sites[i]),
   n_pair_informative_sites=as.integer(ft$n_pair_informative_sites[i]),
   Hudson_Fst=as.numeric(ft$hudson_Fst_signed[i]),BH_q=as.numeric(ft$BH_q_10_pairs[i]),
   CI_low=low,CI_high=hi,interval_available=is.finite(low)&&is.finite(hi)))
}
robust$interpretation <- ifelse(robust$n_pair_informative_sites<2L,
 'UNRELIABLE: fewer than 2 pair-informative sites',
 ifelse(robust$scenario=='282_strict_18_confounded',
 'SENSITIVITY: heterogeneous collection years/studies',
 'EXPLORATORY: one linked gene; not genome-wide differentiation'))
write_tsv(robust,'Table_S34_Fst_sensitivity_and_interval_availability.tsv')

# Correct the earlier geom_segment warning: draw interval only when estimable.
if (requireNamespace('ggplot2', quietly=TRUE)) {
  library(ggplot2)
  main<-robust[robust$scenario=='2020_four_strict',,drop=FALSE]
  main$contrast<-paste(main$state_a,main$state_b,sep=' vs ')
  main$contrast<-factor(main$contrast,levels=rev(main$contrast))
  plot_obj<-ggplot(main,aes(x=Hudson_Fst,y=contrast))+
   geom_vline(xintercept=0,linetype='dashed',color='grey50')+
   geom_segment(data=main[main$interval_available,,drop=FALSE],
     aes(x=CI_low,xend=CI_high,y=contrast,yend=contrast),
     linewidth=.6,color='#426789')+
   geom_point(size=2,color='#B3593E')+
   geom_text(data=main[!main$interval_available,,drop=FALSE],
     aes(label='CI unavailable'),hjust=-0.1,vjust=-.5,size=3.1,color='#794C36')+
   labs(title='PI4K differentiation in the shared-study 2020 cohort',
    subtitle='Four selected linked SNPs; missing bootstrap intervals labelled explicitly',
    x='Hudson-type FST (signed)',y=NULL)+clean_theme()+
   ggplot2::coord_cartesian(clip='off')+
   ggplot2::theme(plot.margin=ggplot2::margin(12,85,12,10))
  export_plot(plot_obj,'Figure_7a_CORRECTED_primary2020_Fst_with_CI_availability',10,6.3)

  pn<-pat_state[pat_state$n>0,,drop=FALSE]
  pn$pattern<-factor(pn$pattern,levels=sort(all_patterns))
  gp<-ggplot(pn,aes(x=state,y=n,fill=pattern))+
     geom_col(width=.72,color='white',linewidth=.2)+
     labs(title='PI4K four-marker consensus patterns in the 2020 cohort',
        subtitle='Complete cases only; patterns are not independently phased parasite haplotypes',
        x='State',y='Samples with all four SNPs called',fill='Allele-code pattern')+
     clean_theme()
  export_plot(gp,'Figure_13_primary2020_four_marker_consensus_patterns',9,5.5)

  lj<-ld[ld$subset=='All_2020',,drop=FALSE]
  # Symmetric site-by-site LD heatmap.
  lp<-rbind(data.frame(site_a=lj$marker_a,site_b=lj$marker_b,r2=lj$r_squared),
            data.frame(site_a=lj$marker_b,site_b=lj$marker_a,r2=lj$r_squared),
            data.frame(site_a=MARKERS,site_b=MARKERS,r2=1))
  lp$site_a<-factor(lp$site_a,levels=MARKERS)
  lp$site_b<-factor(lp$site_b,levels=rev(MARKERS))
  p_ld<-ggplot(lp,aes(site_a,site_b,fill=r2))+geom_tile(color='white')+
       geom_text(aes(label=ifelse(is.finite(r2),sprintf('%.2f',r2),'NA')),size=3.7)+
       scale_fill_gradient(low='white',high='#365D80',limits=c(0,1),na.value='grey85')+
       labs(title='Pairwise LD among four PI4K candidate markers',
         subtitle='r-squared calculated across the 119 consensus-called isolates',
         x='Chromosome 5 SNP position',y='Chromosome 5 SNP position',fill='r-squared')+
       coord_equal()+clean_theme()+theme(panel.grid=element_blank())
  export_plot(p_ld,'Figure_14_primary2020_four_marker_pairwise_LD',7,6)

  cl<-ggplot(call_df,aes(x=state,y=complete_call_fraction))+
   geom_col(fill='#426789',width=.65)+
   geom_errorbar(aes(ymin=wilson_95_low,ymax=wilson_95_high),width=.18,color='#303D47')+
   scale_y_continuous(limits=c(0,1),labels=function(x)paste0(round(100*x),'%'))+
   labs(title='All-four-marker callability by state',
    subtitle='Wilson 95% intervals; only fully called patterns included in pattern counts',
    x='State',y='Fraction of samples with all 4 markers callable')+clean_theme()
  export_plot(cl,'Figure_15_primary2020_four_marker_complete_callability',8,5)
  cat('Exported four figures in TIFF (600 dpi) and PDF.\n')
} else cat('ggplot2 unavailable; numerical tables complete, figures skipped.\n')
cat('STEP 26 COMPLETE: no alterations to input genotype calls or prespecified QC criteria.\n')
