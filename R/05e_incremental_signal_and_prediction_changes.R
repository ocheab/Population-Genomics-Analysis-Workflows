#!/usr/bin/env Rscript
# PI4K Step 30: Descriptive, fixed-prediction audit of incremental genetic information.
# Requires existing Step29A (out-of-fold predictions) and Step29C (selected weights).
# Does NOT refit models, change genomic variants, or claim external validation.
# All tables and scientific figures are generated in R.
options(stringsAsFactors = FALSE)
source('R/01_helpers.R')
cat('STEP 30: Fold-aligned incremental-signal audit (NO model refitting)\n')
D <- file.path(OUT, 'model_step29_corrected')
O <- file.path(OUT, 'model_step30_incremental_audit')
F <- file.path(OUT, 'figures_step30_incremental_audit')
dir.create(O, recursive=TRUE, showWarnings=FALSE)
dir.create(F, recursive=TRUE, showWarnings=FALSE)
read_it <- function(fn) read.delim(required(file.path(D, fn)), sep='\t', quote='"', check.names=FALSE, stringsAsFactors=FALSE)
write_it <- function(x, name) {
  target <- file.path(O, name)
  write.table(x, target, sep='\t', quote=TRUE, row.names=FALSE, na='NA')
  cat('Saved:', target, '\n')
}
P <- read_it('Table_29A_fold_aligned_out_of_fold_predictions.tsv')
H <- read_it('Table_29C_nested_hyperparameters_by_outer_fold.tsv')
M <- read_it('Table_29B_metrics_by_repeat_and_model.tsv')
needed_p <- c('repeat_id','sample_id','truth','outer_fold','QC_only','QC_plus_SNP','SNP_nested')
needed_h <- c('repeat_id','fold','combo_weight_SNP','k_combo_available','n_SNP_abstained','n_test_unseen_alleles')
assert(all(needed_p %in% names(P)), 'Step 29 prediction table is missing required fields')
assert(all(needed_h %in% names(H)), 'Step 29 fold table is missing required fields')
assert(nrow(P)==595L && nrow(H)==25L, 'Expect 119 samples x 5 repeats and 25 outer folds')
assert(!anyDuplicated(P[c('repeat_id','sample_id')]), 'Duplicate predictions')
assert(!anyDuplicated(H[c('repeat_id','fold')]), 'Duplicate hyperparameter fold rows')
assert(all(!is.na(P$QC_only) & !is.na(P$QC_plus_SNP)), 'Unexpected QC/combined abstention')
assert(all(P$truth %in% STATES & P$QC_only %in% STATES & P$QC_plus_SNP %in% STATES),
       'Unexpected state labels in predictions')
reps <- sort(unique(P$repeat_id))
assert(identical(as.integer(reps), 1:5), 'Expect repeats 1 through 5')
ids <- P$sample_id[P$repeat_id==reps[1]]
assert(length(ids)==119L && length(unique(ids))==119L, 'Expected 119 unique samples in first repeat')
for (r in reps) assert(setequal(P$sample_id[P$repeat_id==r], ids), 'Samples differ across repeats')
key <- paste(P$repeat_id, P$outer_fold, sep=':')
kh <- paste(H$repeat_id, H$fold, sep=':')
kmatch <- match(key, kh)
assert(!anyNA(kmatch), 'Could not map all held-out samples to fold weights')
P$snp_weight <- H$combo_weight_SNP[kmatch]
P$k_combo <- H$k_combo_available[kmatch]
P$qc_correct <- P$QC_only==P$truth
P$combined_correct <- P$QC_plus_SNP==P$truth
P$snp_changes_prediction <- P$QC_only!=P$QC_plus_SNP
P$fixes_qc_error <- !P$qc_correct & P$combined_correct
P$breaks_qc_correct <- P$qc_correct & !P$combined_correct
P$both_correct <- P$qc_correct & P$combined_correct
P$both_wrong <- !P$qc_correct & !P$combined_correct
# When weight == 0, scaling the QC scores must preserve the original QC ranking.
assert(all(!P$snp_changes_prediction[P$snp_weight==0]),
       'A zero-SNP-weight fold differs from QC predictions; inspect scoring')
# Pairwise comparison of predictions on exactly the same sample/fold.
by_state_rep <- do.call(rbind, lapply(reps, function(r) {
  z <- P[P$repeat_id==r,,drop=FALSE]
  do.call(rbind,lapply(STATES, function(s) {
    w <- z[z$truth==s,,drop=FALSE]
    n <- nrow(w)
    data.frame(repeat_id=r, state=s, n=n,
               qc_correct=sum(w$qc_correct), combo_correct=sum(w$combined_correct),
               qc_recall=mean(w$qc_correct), combo_recall=mean(w$combined_correct),
               delta_recall=mean(w$combined_correct)-mean(w$qc_correct),
               qc_wrong_fixed=sum(w$fixes_qc_error), qc_correct_broken=sum(w$breaks_qc_correct),
               predictions_changed=sum(w$snp_changes_prediction),
               zero_snp_weight_n=sum(w$snp_weight==0))
  }))
}))
write_it(by_state_rep, 'Table_30A_state_by_repeat_paired_prediction_audit.tsv')
summary_states <- do.call(rbind,lapply(STATES,function(s) {
  t <- by_state_rep[by_state_rep$state==s,,drop=FALSE]
  data.frame(state=s,n_unique_isolates=t$n[1],
             mean_QC_recall=mean(t$qc_recall),mean_combined_recall=mean(t$combo_recall),
             mean_delta_recall=mean(t$delta_recall),
             qc_errors_fixed_over_five_repeats=sum(t$qc_wrong_fixed),
             qc_correct_broken_over_five_repeats=sum(t$qc_correct_broken),
             predictions_changed_over_five_repeats=sum(t$predictions_changed),
             unique_sample_reuse_note='Counts include 5 predictions per isolate; not independent observations')
}))
write_it(summary_states,'Table_30B_state_incremental_recall_and_error_changes.tsv')
# Audit selected SNP weight, number of markers, unseen codes and decision changes by fold.
fold_stats <- do.call(rbind, lapply(seq_len(nrow(H)),function(i) {
  h <- H[i,,drop=FALSE]
  w <- P[P$repeat_id==h$repeat_id & P$outer_fold==h$fold,,drop=FALSE]
  assert(nrow(w)>0L, 'Empty outer fold')
  data.frame(repeat_id=h$repeat_id,fold=h$fold,n_test=nrow(w),
             k_combo_available=h$k_combo_available,
             combo_weight_SNP=h$combo_weight_SNP,
             n_SNP_nested_abstained=h$n_SNP_abstained,
             n_test_unseen_allele_code_occurrences=h$n_test_unseen_alleles,
             note_unseen='Unseen counts from SNP-only panel, not the combined panel',
             n_combined_predictions_changed_vs_QC=sum(w$snp_changes_prediction),
             n_QC_errors_fixed=sum(w$fixes_qc_error),n_QC_correct_broken=sum(w$breaks_qc_correct),
             delta_overall_accuracy=mean(w$combined_correct)-mean(w$qc_correct))
}))
write_it(fold_stats,'Table_30C_fold_level_genetic_weight_and_prediction_effect.tsv')
# Every sample appears in five repeats. Aggregate by sample (not by 595 pseudo-independent rows).
by_sample <- do.call(rbind, lapply(ids, function(id) {
  w <- P[P$sample_id==id,,drop=FALSE]
  data.frame(sample_id=id,state=w$truth[1],n_repeats=nrow(w),
             qc_correct_repeats=sum(w$qc_correct),combo_correct_repeats=sum(w$combined_correct),
             improvement_repeats=sum(w$fixes_qc_error), harm_repeats=sum(w$breaks_qc_correct),
             changed_predictions_repeats=sum(w$snp_changes_prediction),
             folds_weight0=sum(w$snp_weight==0),
             mean_accuracy_delta=mean(w$combined_correct)-mean(w$qc_correct))
}))
write_it(by_sample,'Table_30D_sample_level_prediction_change_stability.tsv')
# Calculate balanced accuracy from state-specific recalls, counting abstentions as errors.
score_BA <- function(truth, pred) {
  mean(vapply(STATES,function(s) {
    w <- which(truth==s)
    if(!length(w)) stop('State absent from resampled evaluation')
    mean(!is.na(pred[w]) & pred[w]==s)
  },numeric(1)))
}
# Predictions matrix: rows are isolates, columns are CV repeats; no fitted model is rerun.
QC <- CO <- matrix(NA_character_,nrow=length(ids),ncol=length(reps))
truth <- P$truth[match(ids,P$sample_id)]
for (j in seq_along(reps)) {
  w <- P[P$repeat_id==reps[j],,drop=FALSE]
  ord <- match(ids,w$sample_id)
  assert(!anyNA(ord) && all(w$truth[ord]==truth), 'Mismatched sample truth or order')
  QC[,j] <- w$QC_only[ord]
  CO[,j] <- w$QC_plus_SNP[ord]
}
mean_delta <- function(ix) {
  mean(vapply(seq_along(reps),function(j)
    score_BA(truth[ix],CO[ix,j])-score_BA(truth[ix],QC[ix,j]),numeric(1)))
}
observed <- mean_delta(seq_along(ids))
reported <- mean(subset(M,model=='QC_plus_SNP')$balanced_accuracy) -
            mean(subset(M,model=='QC_only')$balanced_accuracy)
assert(abs(observed-reported)<1e-10, 'Matched-prediction BA delta differs from Step 29 reported value')
# Patient/isolate-level stratified bootstrap: the SAME sampled isolate indices are
# used for all five repeats and both models, preserving repeated-measure dependence.
# The fits/predictions remain fixed: this is NOT a refit-bootstrap CI and is not
# an inferential test of conditional genomic information.
B <- as.integer(Sys.getenv('PI4K_STEP30_BOOT','2000'))
assert(is.finite(B) && B>=100L && B<=20000L,'PI4K_STEP30_BOOT must be between 100 and 20000')
set.seed(20261109L)
state_index <- lapply(STATES,function(s) which(truth==s))
boot <- numeric(B)
for(b in seq_len(B)) {
  ix <- unlist(lapply(state_index,function(w) sample(w,size=length(w),replace=TRUE)),use.names=FALSE)
  boot[b] <- mean_delta(ix)
}
boot_ci <- unname(quantile(boot,c(.025,.5,.975),names=FALSE,na.rm=TRUE))
write_it(data.frame(replicate=seq_len(B),delta_balanced_accuracy=boot),
         'Table_30E_fixed_prediction_isolate_bootstrap_distribution.tsv')
write_it(data.frame(n_unique_isolates=length(ids),n_repeats=length(reps),n_bootstrap=B,
                    mean_delta_balanced_accuracy=observed,
                    percentile_2_5=boot_ci[1],bootstrap_median=boot_ci[2],percentile_97_5=boot_ci[3],
                    interpretation='Conditional on saved OOF predictions: no model refits; NOT a definitive generalisation CI or SNP-adjusted-value hypothesis test'),
         'Table_30F_fixed_prediction_bootstrap_summary.tsv')
# Leave-one-isolate-out influence across all five repeats; no refitting.
loo <- do.call(rbind,lapply(seq_along(ids),function(i) {
  val <- mean_delta(setdiff(seq_along(ids),i))
  data.frame(sample_id=ids[i],state=truth[i],
             delta_if_omitted=val,impact_on_delta=observed-val)
}))
write_it(loo[order(-abs(loo$impact_on_delta)),,drop=FALSE],
         'Table_30G_isolate_influence_no_refit.tsv')
checks <- data.frame(check=c('119 unique held-out isolates',
                              'Five correlated repeats cover identical isolates',
                              'Every held-out row maps to fold SNP weight',
                              'Lambda-zero folds reproduce QC-only predictions',
                              'No abstentions in QC or combined predictions',
                              'Reported Step29 balanced-accuracy delta reproduced',
                              'Allele-code events are from SNP-only panel, not confirmed combined-panel fallback'),
                     status='PASS',stringsAsFactors=FALSE)
write_it(checks,'Table_30H_internal_consistency_checks.tsv')
writeLines(c(
  'STEP 30 — Interpretation safeguards',
  'Data: 119 isolates, 5 repeated cross-validation predictions, 25 outer test folds.',
  'Do not treat 595 held-out predictions as 595 independent individuals.',
  'The stratified bootstrap resamples 119 isolates, retaining their five predictions.',
  'All fitted models and held-out predictions remain FIXED during bootstrap.',
  'Result is a conditional descriptive percentile interval, NOT a refit-bootstrap CI.',
  'No p-value for incremental SNP information conditional on QC is calculated.',
  'No combined-model fallback flag was saved by Step29.',
  'A combined prediction identical to QC DOES NOT prove fallback.',
  'SNP-only abstentions/unseen-code events are not necessarily the combined-panel events.',
  'At SNP weight 0, combined and QC predicted labels must agree.',
  'Five-state assignment is exploratory and has no independent geographic validation.',
  sprintf('Observed combined minus QC mean balanced-accuracy difference: %.8f',observed)),
  file.path(O,'Step30_interpretation_limitations.txt'))
if(requireNamespace('ggplot2',quietly=TRUE)) {
  library(ggplot2)
  save_fig <- function(plot,name,width=9,height=5.5) {
    ggplot2::ggsave(file.path(F,paste0(name,'.pdf')),plot=plot,width=width,height=height,
                    units='in',device='pdf',useDingbats=FALSE)
    ggplot2::ggsave(file.path(F,paste0(name,'.tiff')),plot=plot,width=width,height=height,
                    units='in',dpi=600,compression='lzw')
  }
  state_long <- rbind(data.frame(state=summary_states$state,model='QC only',
                                recall=summary_states$mean_QC_recall),
                      data.frame(state=summary_states$state,model='QC plus PI4K',
                                 recall=summary_states$mean_combined_recall))
  p1 <- ggplot(state_long,aes(x=state,y=recall,group=model,color=model))+
    geom_point(size=3)+geom_line(linewidth=.8)+
    scale_y_continuous(limits=c(0,1),labels=function(x)paste0(round(x*100),'%'))+
    labs(title='State-level recall from matched held-out predictions',
         subtitle='Means across five repeated splits of the same 119 isolates',
         x='Collection state',y='Recall',color='Model')+clean_theme()
  save_fig(p1,'Figure_30A_state_recall_QC_vs_combined')
  change_long <- rbind(data.frame(state=summary_states$state,
                                 outcome='QC errors corrected',n=summary_states$qc_errors_fixed_over_five_repeats),
                       data.frame(state=summary_states$state,
                                  outcome='QC correct made incorrect',n=summary_states$qc_correct_broken_over_five_repeats))
  p2 <- ggplot(change_long,aes(x=state,y=n,fill=outcome))+
    geom_col(position='dodge')+
    labs(title='How PI4K changed QC-only decisions',
         subtitle='Counts across five repeated held-out predictions per isolate (not independent)',
         x='Collection state',y='Prediction instances',fill='Change')+clean_theme()
  save_fig(p2,'Figure_30B_QC_prediction_corrections_and_harms')
  p3 <- ggplot(data.frame(delta=boot),aes(x=100*delta))+
    geom_histogram(bins=35,fill='grey55',color='white')+
    geom_vline(xintercept=100*observed,linetype='dashed')+
    labs(title='Isolate-resampled sensitivity of observed model difference',
         subtitle='Models held fixed; conditional descriptive bootstrap, not external-validation CI',
         x='QC+PI4K minus QC balanced accuracy (percentage points)',y='Bootstrap replicates')+clean_theme()
  save_fig(p3,'Figure_30C_fixed_prediction_isolate_bootstrap')
  cat('Figures exported: 3 PDF and 3 600-dpi TIFF files.\n')
} else cat('ggplot2 unavailable; numeric tables completed, figures skipped.\n')
cat(sprintf('STEP 30 COMPLETE. Combined minus QC BA: +%.4f (%.2f percentage points)\n',
            observed,100*observed))
cat(sprintf('Conditional fixed-prediction percentile interval: %.2f to %.2f percentage points\n',
            100*boot_ci[1],100*boot_ci[3]))
cat('CAVEAT: This is not an independently validated barcode or formal conditional-SNP test.\n')
