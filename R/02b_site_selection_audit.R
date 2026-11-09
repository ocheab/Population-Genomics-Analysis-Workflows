# Run from PI4K_Southwest_R_Analysis_Suite directory, after 02_preflight_qc.R.
# No source VCF modification or genomic download.
source('R/01_helpers.R')
cat('WP2b: audited site-level filtering and state-wise allele counts\n')
cohorts <- readRDS(required(file.path(OUT, 'rds', 'QC_cohorts.rds')))

summary_rows <- list()
marker_rows <- list()
state_rows <- list()
for (cohort in c('primary2020','all282')) {
  z <- cohorts[[cohort]]
  G <- z$X$G
  ss <- z$site
  m <- z$meta
  assert(ncol(G)==nrow(m) && nrow(G)==nrow(ss),'cohort dimensions disagree')
  assert(identical(colnames(G), m$sample_id),'Sample-ID order mismatch')
  assert(identical(rownames(G), ss$id),'Site row order mismatch')
  cand <- which(ss$allele_count==2L)
  mark <- which(ss$informative)
  cat(sprintf('%s: %s consensus polymorphic, %s prespecified high-quality markers\n',cohort,length(cand),length(mark)))
  # Each threshold is independently evaluated. Do not portray these as a
  # sequential attrition cascade, because the filters overlap.
  for (crit in c('all_PASS','consensus_polymorphic','biallelic_MAC_ge3',
                 'overall_call_ge90','each_state_call_ge80','all_prespecified')) {
    keep <- switch(crit,
      all_PASS=rep(TRUE,nrow(ss)),
      consensus_polymorphic=ss$allele_count==2L,
      biallelic_MAC_ge3=ss$allele_count==2L & !is.na(ss$minor_count) & ss$minor_count>=3,
      overall_call_ge90=ss$callable>=SETTINGS$call_rate,
      each_state_call_ge80=rowSums(ss[,STATES,drop=FALSE]>=SETTINGS$state_call_rate)==length(STATES),
      all_prespecified=ss$informative)
    summary_rows[[length(summary_rows)+1L]] <- data.frame(cohort=cohort,
      criterion=crit, sites=sum(keep,na.rm=TRUE))
  }
  for (i in mark) {
    marker_rows[[length(marker_rows)+1L]] <- data.frame(cohort=cohort,
      chrom=z$X$sites$chrom[i],pos_1based=ss$pos[i], ref=ss$ref[i],alt=ss$alt[i],
      n_consensus_alleles=ss$allele_count[i],minor_consensus_sample_count=ss$minor_count[i],
      overall_call_fraction=ss$callable[i],
      Oyo_call_fraction=ss$Oyo[i],Ogun_call_fraction=ss$Ogun[i],
      Lagos_call_fraction=ss$Lagos[i],Osun_call_fraction=ss$Osun[i],
      Ondo_call_fraction=ss$Ondo[i])
  }
  # Within-state fixed sites (all REF or all ALT) must be retained as rows.
  for(i in cand) {
    global_codes <- sort(unique(as.integer(na.omit(G[i,]))))
    if(!identical(global_codes,c(0L,1L))) next # describe multiallelic separately
    for(s in STATES) {
      x <- G[i,m$state==s]
      obs <- x[!is.na(x)]
      n0 <- sum(obs==0L); n1 <- sum(obs==1L)
      assert(n0+n1==length(obs),'Unexpected allele code at REF/ALT1 site')
      state_rows[[length(state_rows)+1L]] <- data.frame(cohort=cohort,
        chrom=z$X$sites$chrom[i],pos_1based=ss$pos[i],state=s,
        n_samples=sum(m$state==s),n_consensus_called=length(obs),
        n_missing=sum(is.na(x)),n_REF=n0,n_ALT1=n1,
        ALT1_frequency=if(length(obs)>0) n1/length(obs) else NA_real_,
        prespecified_high_quality=ss$informative[i])
    }
  }
}
write_tsv(do.call(rbind,summary_rows),'Table_S19_site_QC_overlapping_filter_counts.tsv')
if(length(marker_rows)) write_tsv(do.call(rbind,marker_rows),
  'Table_S20_prespecified_high_quality_SNP_positions.tsv')
if(length(state_rows)) write_tsv(do.call(rbind,state_rows),
  'Table_S21_state_REF_ALT1_counts_all_consensus_biallelic_sites.tsv')
cat('WP2b COMPLETE. Reviewed marker counts, overlapping site-QC filters, and within-state fixed calls.\n')
