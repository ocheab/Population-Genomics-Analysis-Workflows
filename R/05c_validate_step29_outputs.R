#!/usr/bin/env Rscript
# Step 29 post-run numerical and data-isolation checks (no genotype re-fitting).
source('R/01_helpers.R')
cohort <- readRDS(required(file.path(OUT,'rds','QC_cohorts.rds')))$primary2020
D <- file.path(OUT,'model_step29_corrected')
r <- function(fn) read.delim(required(file.path(D,fn)),check.names=FALSE,
                              stringsAsFactors=FALSE,quote='"')
P <- r('Table_29A_fold_aligned_out_of_fold_predictions.tsv')
M <- r('Table_29B_metrics_by_repeat_and_model.tsv')
H <- r('Table_29C_nested_hyperparameters_by_outer_fold.tsv')
models <- c('majority','QC_only','SNP_nested','QC_plus_SNP','SNP_fixed4_posthoc')
assert(all(models %in% names(P)), 'Missing model predictions')
assert(nrow(P)==119L*5L, 'Expected five repeats of 119 samples')
assert(!anyDuplicated(P[c('repeat_id','sample_id')]), 'Duplicate sample/repeat rows')
assert(nrow(M)==5L*length(models), 'Unexpected metric rows')
assert(nrow(H)==25L, 'Expected 25 outer-fold hyperparameter rows')
meta <- cohort$meta
assert(setequal(unique(P$sample_id),meta$sample_id), 'Prediction/sample ID mismatch')
original_case <- as.character(meta$case_id)
original_case[is.na(original_case)|original_case==''] <-
  as.character(meta$sample_id[is.na(original_case)|original_case==''])
for(rep in sort(unique(P$repeat_id))) {
  z <- P[P$repeat_id==rep,,drop=FALSE]
  assert(nrow(z)==119L && setequal(z$sample_id,meta$sample_id),'Incomplete repeat')
  iy <- match(z$sample_id,meta$sample_id)
  assert(identical(z$truth,as.character(meta$state[iy])), 'Truth labels do not match cohort metadata')
  case <- original_case[iy]
  assert(all(vapply(split(z$outer_fold,case),function(f)length(unique(f))==1L,logical(1))),
         'A recorded case group crossed outer folds')
  assert(length(unique(z$outer_fold))==5L,'Not five distinct outer folds')
  for(m in models) {
    pred <- z[[m]]
    assert(all(is.na(pred)|pred %in% STATES),paste('Invalid predicted label',m))
    tp <- vapply(STATES,function(s)sum(z$truth==s & !is.na(pred) & pred==s),numeric(1))
    ns <- vapply(STATES,function(s)sum(z$truth==s),numeric(1))
    rec <- tp/ns
    pre <- vapply(STATES,function(s) {den<-sum(!is.na(pred)&pred==s);if(den)tp[s]/den else 0},numeric(1))
    f1 <- ifelse(pre+rec>0,2*pre*rec/(pre+rec),0)
    obs <- c(balanced_accuracy=mean(rec),macro_F1=mean(f1),
             overall_accuracy=mean(!is.na(pred)&z$truth==pred),coverage=mean(!is.na(pred)))
    row <- M[M$repeat_id==rep & M$model==m,,drop=FALSE]
    assert(nrow(row)==1L,paste('Missing metric row:',m,rep))
    err <- max(abs(as.numeric(unlist(row[names(obs)],use.names=FALSE))-unname(obs)))
    assert(is.finite(err)&&err<1e-10,paste('Metrics disagree for',m,rep,err))
  }
}
assert(all(H$k_snp_requested>=1L & H$k_combo_requested>=1L),'Invalid requested panel size')
assert(all(H$combo_weight_SNP>=0 & H$combo_weight_SNP<=1),'Invalid combined SNP weight')
assert(all(H$k_snp_available<=H$k_snp_requested),'Available SNP count exceeds requested')
assert(all(H$k_combo_available<=H$k_combo_requested),'Available combined SNP count exceeds requested')
report <- data.frame(
 check=c('Predictions unique by sample and repeat','All 119 samples represented in five repeats',
         'Truth labels agree with original 2020 metadata','Case groups isolated within folds',
         'All held-out performance metrics exactly reproduced',
         'All hyperparameters within configured bounds'),
 result='PASS',stringsAsFactors=FALSE)
write.table(report,file.path(D,'Table_29L_post_run_validation_checks.tsv'),
  sep='\t',row.names=FALSE,quote=TRUE)
cat('STEP 29 POST-RUN VALIDATION PASSED — six checks.\n')
