#!/usr/bin/env Rscript
# Step 29 plot-only finalization. Reuses existing held-out results; DOES NOT fit models.
# Run from PI4K_Southwest_R_Analysis_Suite project root, after Step 29 model tables exist.
source('R/01_helpers.R')
result_dir <- file.path(OUT, 'model_step29_corrected')
figure_dir <- file.path(OUT, 'figures_step29_corrected')
dir.create(figure_dir, recursive=TRUE, showWarnings=FALSE)
read_table <- function(name) {
  f <- file.path(result_dir, name)
  stopifnot(file.exists(f))
  read.delim(f, sep='\t', quote='"', header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
}
MET <- read_table('Table_29B_metrics_by_repeat_and_model.tsv')
M <- read_table('Table_29G_fold_aligned_paired_descriptive_deltas.tsv')
state_summary <- read_table('Table_29I_state_recall_and_coverage.tsv')
mean_perf <- read_table('Table_29F_repeated_mean_performance.tsv')
stopifnot(nrow(M) == 5L, length(unique(M$repeat_id)) == 5L)
stopifnot(all(c('repeat_id','model','balanced_accuracy') %in% names(MET)))
stopifnot(all(c('state','model','mean_repeated_recall') %in% names(state_summary)))
stopifnot(all(c('delta_SNP_nested_minus_QC_balanced_accuracy',
                'delta_QC_plus_SNP_minus_QC_balanced_accuracy') %in% names(M)))
stopifnot(requireNamespace('ggplot2', quietly=TRUE))
library(ggplot2)
export_one <- function(p, stem, w=9, h=5) {
  pdf_file <- file.path(figure_dir, paste0(stem, '.pdf'))
  tif_file <- file.path(figure_dir, paste0(stem, '.tiff'))
  ggsave(pdf_file, plot=p, width=w, height=h, units='in', device='pdf', useDingbats=FALSE)
  ggsave(tif_file, plot=p, width=w, height=h, units='in', dpi=600, compression='lzw')
  stopifnot(file.exists(pdf_file), file.exists(tif_file),
            file.info(pdf_file)$size > 0, file.info(tif_file)$size > 0)
  c(pdf_file, tif_file)
}
p <- ggplot(MET, aes(x=factor(repeat_id), y=balanced_accuracy, color=model, group=model)) +
  geom_hline(yintercept=0.2, linetype='dashed', colour='grey65') +
  geom_line(linewidth=.65) + geom_point(size=2) +
  scale_y_continuous(limits=c(0,1), labels=function(x) sprintf('%d%%', round(x*100))) +
  labs(title='PI4K genetic and QC model performance, matched folds',
       subtitle='Five repeated nested CV partitions of the same 119 isolates',
       x='Outer fold repeat', y='Balanced accuracy', color='Model') + clean_theme()
f1 <- export_one(p, 'Figure_29A_matched_fold_balanced_accuracy')
dd <- data.frame(repeat_id=M$repeat_id,
                 snp=M$delta_SNP_nested_minus_QC_balanced_accuracy,
                 combined=M$delta_QC_plus_SNP_minus_QC_balanced_accuracy)
long <- rbind(data.frame(repeat_id=dd$repeat_id, model='SNP minus QC', delta=dd$snp),
              data.frame(repeat_id=dd$repeat_id, model='Combined minus QC', delta=dd$combined))
p2 <- ggplot(long, aes(x=factor(repeat_id), y=100*delta, group=model, color=model)) +
  geom_hline(yintercept=0, linetype='dashed', colour='grey60') +
  geom_point(size=2.3) + geom_line(linewidth=.65) +
  labs(title='Paired descriptive change in held-out balanced accuracy',
       subtitle='Matched folds; not five independent studies',
       x='Outer fold repeat', y='Balanced-accuracy difference (percentage points)',
       color='Comparison') + clean_theme()
f2 <- export_one(p2, 'Figure_29B_model_minus_QC_differences')
p3 <- ggplot(state_summary, aes(x=state, y=mean_repeated_recall, fill=model)) +
  geom_col(position='dodge') +
  scale_y_continuous(limits=c(0,1), labels=function(x) sprintf('%d%%', round(x*100))) +
  labs(title='Class-specific recall, repeated held-out predictions',
       subtitle='Abstentions counted as incorrect; repeated splits reuse isolates',
       x='Collection state', y='Mean recall', fill='Model') + clean_theme()
f3 <- export_one(p3, 'Figure_29C_per_state_recall')
files <- c(f1,f2,f3)
summary_row <- function(m) mean_perf[mean_perf$model==m, , drop=FALSE]
stopifnot(nrow(summary_row('QC_only'))==1L,
          nrow(summary_row('QC_plus_SNP'))==1L)
delta <- summary_row('QC_plus_SNP')$balanced_accuracy -
         summary_row('QC_only')$balanced_accuracy
report <- data.frame(item=c('Plot export','Six PDF/TIFF files','Mean combined-minus-QC balanced accuracy',
                              'Source of values','Interpretation'),
                     value=c('PASS', as.character(length(files)),
                             sprintf('%.9f',delta),
                             'Existing Step 29 out-of-fold summary tables',
                             'Exploratory; correlated repeated CV; no external validation'))
write.table(report, file.path(result_dir,'Table_29N_figure_completion_audit.tsv'),
            sep='\t', row.names=FALSE, quote=TRUE)
cat('STEP 29 FIGURE-ONLY FINALIZATION COMPLETE: three PDFs and three 600-dpi TIFFs.\n')
cat(sprintf('Descriptive QC+SNP minus QC balanced accuracy: %+0.4f percentage points\n',100*delta))
cat('No models refit; no input genotype or QC thresholds changed.\n')
