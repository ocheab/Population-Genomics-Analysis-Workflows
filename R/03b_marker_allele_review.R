# STEP 24: Verify allele identity, state-level counts and figures at candidate PI4K SNPs.
# Run from PI4K_Southwest_R_Analysis_Suite after R/02_preflight_qc.R.
# Does not read remote data or change source VCFs or existing RDS.
source('R/01_helpers.R')
cat('STEP 24: PI4K allele-code and state-frequency audit (local RDS only)\n')
cohorts <- readRDS(required(file.path(OUT, 'rds', 'QC_cohorts.rds')))
stopifnot(all(c('primary2020', 'all282') %in% names(cohorts)))

site_info <- list(); state_info <- list()
for (cohort in c('primary2020', 'all282')) {
  z <- cohorts[[cohort]]
  G <- z$X$G; ss <- z$site; m <- z$meta
  assert(nrow(G)==nrow(ss) && ncol(G)==nrow(m), 'Matrix dimensions mismatch')
  assert(identical(colnames(G), as.character(m$sample_id)), 'Sample order differs from QC metadata')
  assert(identical(rownames(G), as.character(ss$id)), 'Site order differs from QC metadata')
  assert(all(as.character(m$state) %in% STATES), 'Unexpected state code')
  ix <- which(ss$allele_count==2L)
  for (i in ix) {
    called <- as.integer(G[i, !is.na(G[i,])])
    codes <- sort(unique(called))
    assert(length(codes)==2L, paste('Unexpected allele count at',ss$pos[i]))
    allele_names <- c(as.character(ss$ref[i]), strsplit(as.character(ss$alt[i]), ',', fixed=TRUE)[[1]])
    assert(all(codes>=0L & codes < length(allele_names)),
      paste('Allele code outside VCF REF/ALT enumeration at',ss$pos[i]))
    counts <- vapply(codes, function(x) sum(called==x), integer(1))
    # Deterministic tie handling: lowest VCF code is the nominal minor allele when tied.
    minor_code <- codes[which.min(counts)]
    other_code <- codes[which.max(counts)]
    is_selected <- isTRUE(ss$informative[i])
    ref_alt_multiple <- grepl(',', as.character(ss$alt[i]), fixed=TRUE)
    site_info[[length(site_info)+1L]] <- data.frame(
      cohort=cohort, chrom=as.character(z$X$sites$chrom[i]),
      pos_1based=as.integer(ss$pos[i]), REF=as.character(ss$ref[i]),
      VCF_ALT=as.character(ss$alt[i]), multi_ALT_record=ref_alt_multiple,
      observed_code_1=codes[1], observed_allele_1=allele_names[codes[1]+1L],
      observed_code_2=codes[2], observed_allele_2=allele_names[codes[2]+1L],
      minor_allele_code=minor_code, minor_allele=allele_names[minor_code+1L],
      pooled_minor_n=sum(called==minor_code), pooled_n_called=length(called),
      pooled_minor_fraction=sum(called==minor_code)/length(called),
      overall_call_fraction=mean(!is.na(G[i,])),
      high_quality_all_prespecified=is_selected,
      stringsAsFactors=FALSE)
    for (s in STATES) {
      x <- as.integer(G[i, m$state==s])
      n_total <- length(x); n_called <- sum(!is.na(x))
      n_minor <- sum(x==minor_code,na.rm=TRUE)
      freq <- if(n_called) n_minor/n_called else NA_real_
      ci <- wilson_ci(n_minor,n_called)
      state_info[[length(state_info)+1L]] <- data.frame(
        cohort=cohort, state=s, pos_1based=as.integer(ss$pos[i]),
        REF=as.character(ss$ref[i]), VCF_ALT=as.character(ss$alt[i]),
        observed_allele_1=allele_names[codes[1]+1L],
        observed_allele_2=allele_names[codes[2]+1L],
        pooled_minor_allele=allele_names[minor_code+1L],
        pooled_minor_code=minor_code, n_state=n_total, n_called=n_called,
        n_missing=n_total-n_called, call_fraction=n_called/n_total,
        n_pooled_minor=n_minor, pooled_minor_frequency=freq,
        wilson_low=ci[1], wilson_high=ci[2],
        high_quality_all_prespecified=is_selected,
        stringsAsFactors=FALSE)
    }
  }
}
SITES <- do.call(rbind,site_info)
FREQ <- do.call(rbind,state_info)
assert(nrow(SITES[SITES$cohort=='primary2020',])==43L,
  'Expected 43 consensus-polymorphic 2020 sites; please investigate rather than assume')
assert(nrow(SITES[SITES$cohort=='all282',])==77L,
  'Expected 77 consensus-polymorphic 282-sample sites; please investigate rather than assume')
assert(sum(SITES$cohort=='primary2020' & SITES$high_quality_all_prespecified)==4L,
  'Expected four 2020 markers')
assert(sum(SITES$cohort=='all282' & SITES$high_quality_all_prespecified)==18L,
  'Expected 18 all282 markers')
assert(all(FREQ$n_called+FREQ$n_missing==FREQ$n_state),'State call accounting mismatch')
assert(all(is.na(FREQ$pooled_minor_frequency)|
  (FREQ$pooled_minor_frequency>=0 & FREQ$pooled_minor_frequency<=1)), 'Frequency outside 0..1')
# Important: each site's minor allele is defined in the matching cohort, so
# allele-frequency comparisons across cohorts should use named nucleotides.
p1 <- write_tsv(SITES,'Table_S22_all_consensus_polymorphic_allele_code_audit.tsv')
p2 <- write_tsv(FREQ,'Table_S23_state_named_allele_frequencies_and_Wilson_CI.tsv')
p3 <- write_tsv(SITES[SITES$high_quality_all_prespecified,,drop=FALSE],
  'Table_S24_high_quality_marker_nucleotide_identity.tsv')
p4 <- write_tsv(FREQ[FREQ$high_quality_all_prespecified,,drop=FALSE],
  'Table_S25_high_quality_marker_state_frequencies.tsv')
cat('Saved:',basename(p1),basename(p2),basename(p3),basename(p4),'\n')
cat('Multi-ALT VCF records among candidate markers:',
    sum(SITES$multi_ALT_record & SITES$high_quality_all_prespecified), '\n')
cat('High-quality markers with observed code other than 0/1:',
    sum(SITES$high_quality_all_prespecified &
          !(SITES$observed_code_1==0L & SITES$observed_code_2==1L)), '\n')

if(!requireNamespace('ggplot2',quietly=TRUE)) {
  cat('FIGURES NOT EXPORTED: install ggplot2 with:\n')
  cat('Rscript -e \'install.packages("ggplot2", repos="https://cloud.r-project.org")\'\n')
} else {
  library(ggplot2)
  make_allele_plot <- function(which_cohort, number, height, columns) {
    sub <- FREQ[FREQ$cohort==which_cohort & FREQ$high_quality_all_prespecified,,drop=FALSE]
    assert(nrow(sub)>0,'No high-quality markers for plotting')
    # Distinguish alternate nucleotide if the VCF ALT field is multi-allelic.
    sub$marker <- paste0(sub$pos_1based,' (',sub$pooled_minor_allele,')')
    sub$state <- factor(sub$state, levels=STATES)
    sub$marker <- factor(sub$marker, levels=unique(sub$marker))
    plot <- ggplot(sub,aes(x=state,y=pooled_minor_frequency,color=state))+
      geom_errorbar(aes(ymin=wilson_low,ymax=wilson_high),width=.16,alpha=.8)+
      geom_point(size=2.2)+
      facet_wrap(~marker,ncol=columns)+
      scale_y_continuous(limits=c(0,1),labels=function(x)paste0(round(x*100),'%'))+
      labs(title=paste0('PI4K candidate-SNP frequencies — ',which_cohort),
        subtitle='Cohort-defined minor nucleotide; Wilson 95% intervals; descriptive only',
        x=NULL,y='Minor nucleotide frequency')+
      clean_theme()+theme(axis.text.x=element_text(angle=40,hjust=1),legend.position='none')
    export_plot(plot,number,w=if(which_cohort=='primary2020')9 else 12,h=height)
  }
  make_allele_plot('primary2020','Figure_10_primary2020_four_marker_frequencies',6.5,2)
  make_allele_plot('all282','Figure_11_all282_18_marker_frequencies',13,3)
  # 2020 marker-specific missingness, explicitly showing Oyo callability.
  sub <- FREQ[FREQ$cohort=='primary2020' & FREQ$high_quality_all_prespecified,,drop=FALSE]
  sub$state <- factor(sub$state,levels=STATES)
  sub$marker <- factor(as.character(sub$pos_1based),levels=unique(as.character(sub$pos_1based)))
  p <- ggplot(sub,aes(x=state,y=call_fraction,fill=state))+
    geom_col(width=.7)+facet_wrap(~marker,ncol=2)+
    geom_hline(yintercept=0.8,linetype='dashed',color='gray40')+
    scale_y_continuous(limits=c(0,1),labels=function(x)paste0(round(x*100),'%'))+
    labs(title='Call rate at four 2020 PI4K candidate SNPs',
      subtitle='Dashed line: 80% per-state minimum; markers are not validated barcodes',
      x=NULL,y='Call rate')+clean_theme()+theme(legend.position='none',axis.text.x=element_text(angle=40,hjust=1))
  export_plot(p,'Figure_12_primary2020_marker_callability_by_state',9,6.5)
  cat('Three R figures exported as 600-dpi LZW TIFF and vector PDF.\n')
}
cat('STEP 24 COMPLETE. NO input genotypes or thresholds altered.\n')
