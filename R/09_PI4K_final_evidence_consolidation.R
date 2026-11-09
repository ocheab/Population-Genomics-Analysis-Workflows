#!/usr/bin/env Rscript
# Step 31: publication-oriented consolidation of the already completed PI4K analyses.
# NO model fitting, NO genotype alteration, NO new statistical significance test.
# Run in PI4K_Southwest_R_Analysis_Suite, with project directory one level above.
options(stringsAsFactors = FALSE)
project <- normalizePath('..', mustWork=TRUE)
base <- file.path(project, 'results', 'Southwest5', 'R_Analysis_v1')
dest <- file.path(base, 'Step31_publication_evidence')
fig <- file.path(dest, 'figures')
dir.create(dest, showWarnings=FALSE, recursive=TRUE)
dir.create(fig, showWarnings=FALSE, recursive=TRUE)

get_tab <- function(relative) {
  path <- file.path(base, relative)
  if (!file.exists(path)) stop('Required Step 31 input missing: ', path)
  x <- read.delim(path, header=TRUE, sep='\t', quote='"', check.names=FALSE,
                  stringsAsFactors=FALSE, na.strings=c('NA',''))
  attr(x, 'source_file') <- path
  x
}
require_cols <- function(d, keys) {
  missing <- setdiff(keys, names(d))
  if(length(missing)) stop('Missing required columns: ', paste(missing,collapse=', '))
}
write_tab <- function(d, name) {
  path <- file.path(dest, name)
  write.table(d, path, sep='\t', quote=TRUE, row.names=FALSE, na='NA')
  cat('Saved: ', path, '\n', sep='')
}

# Previously audited outputs: source records remain unchanged.
met <- get_tab('model_step29_corrected/Table_29F_repeated_mean_performance.tsv')
delta <- get_tab('model_step30_incremental_audit/Table_30F_fixed_prediction_bootstrap_summary.tsv')
state <- get_tab('model_step30_incremental_audit/Table_30B_state_incremental_recall_and_error_changes.tsv')
check <- get_tab('model_step30_incremental_audit/Table_30H_internal_consistency_checks.tsv')
post <- get_tab('model_step29_corrected/Table_29L_post_run_validation_checks.tsv')
fst <- get_tab('tables/Table_S34_Fst_sensitivity_and_interval_availability.tsv')
markers <- get_tab('tables/Table_S28_four_marker_global_state_exact_tests.tsv')
cohort <- get_tab('tables/Table_1_study_population_and_QC.tsv')

require_cols(met, c('model','balanced_accuracy','macro_F1','overall_accuracy','coverage'))
require_cols(delta, c('n_unique_isolates','n_repeats','n_bootstrap',
                      'mean_delta_balanced_accuracy','percentile_2_5','percentile_97_5'))
require_cols(state, c('state','n_unique_isolates','mean_QC_recall','mean_combined_recall',
  'qc_errors_fixed_over_five_repeats','qc_correct_broken_over_five_repeats',
  'predictions_changed_over_five_repeats'))
require_cols(fst, c('state_a','state_b','scenario','n_selected_sites',
  'n_pair_informative_sites','Hudson_Fst','BH_q','interval_available'))
require_cols(markers, c('pos_1based','p_Fisher_montecarlo_5x2','BH_q_four_markers'))

if (nrow(delta)!=1L || delta$n_unique_isolates[1]!=119 || delta$n_repeats[1]!=5) {
  stop('Unexpected isolate or repeated-CV counts; review changed upstream results')
}
if (!all(c('QC_only','QC_plus_SNP','SNP_nested','SNP_fixed4_posthoc','majority') %in% met$model)) {
  stop('Required Step29 comparator absent')
}
if (!all(as.character(check$status)=='PASS') || !all(as.character(post$result)=='PASS')) {
  stop('Upstream internal consistency gate has not passed')
}
getmod <- function(name, col) as.numeric(met[met$model==name, col][1])
qc <- getmod('QC_only','balanced_accuracy')
combo <- getmod('QC_plus_SNP','balanced_accuracy')
snp <- getmod('SNP_nested','balanced_accuracy')
fixed <- getmod('SNP_fixed4_posthoc','balanced_accuracy')
observed <- as.numeric(delta$mean_delta_balanced_accuracy[1])
lo <- as.numeric(delta$percentile_2_5[1]); hi <- as.numeric(delta$percentile_97_5[1])
if (abs(combo-qc-observed)>1e-6) stop('Step29 and Step30 deltas disagree')
if (!identical(sort(as.character(state$state)), sort(c('Oyo','Ogun','Lagos','Osun','Ondo')))) {
  stop('State set unexpected')
}
fixes <- sum(state$qc_errors_fixed_over_five_repeats)
breaks <- sum(state$qc_correct_broken_over_five_repeats)
changes <- sum(state$predictions_changed_over_five_repeats)
net <- fixes-breaks
if (net != 8L || changes!=232L) stop('State-level transition counts changed; re-audit')
expected_accuracy_delta <- net/(119*5)
if (abs((getmod('QC_plus_SNP','overall_accuracy')-getmod('QC_only','overall_accuracy'))-
        expected_accuracy_delta)>1e-6) stop('Overall accuracy net correction check failed')

model_out <- met
model_out$balanced_accuracy_pct <- 100*as.numeric(model_out$balanced_accuracy)
model_out$macro_F1_pct <- 100*as.numeric(model_out$macro_F1)
model_out$overall_accuracy_pct <- 100*as.numeric(model_out$overall_accuracy)
model_out$coverage_pct <- 100*as.numeric(model_out$coverage)
model_out$analysis_status <- ifelse(model_out$model=='SNP_fixed4_posthoc',
   'POST-HOC; not independently validated','Repeated nested internal CV only')
write_tab(model_out, 'Table_31A_model_performance.tsv')

state_out <- state
state_out$mean_delta_recall <- state$mean_combined_recall-state$mean_QC_recall
state_out$net_correct_across_five_repeats <- state$qc_errors_fixed_over_five_repeats -
  state$qc_correct_broken_over_five_repeats
state_out$note <- 'Five correlated predictions per isolate: counts are not independent isolates'
write_tab(state_out, 'Table_31B_geographic_prediction_tradeoffs.tsv')

fst_out <- fst
fst_out$BH_q <- as.numeric(fst_out$BH_q)
fst_out$site_caution <- ifelse(fst_out$n_pair_informative_sites<2,
  'NOT RELIABLY ESTIMABLE', 'Single linked PI4K gene; exploratory')
fst_out$significance_caution <- ifelse(fst_out$BH_q<0.05,
  'BH q<0.05 within scenario; NOT an independently validated barcode',
  'BH q>=0.05 within scenario')
write_tab(fst_out, 'Table_31C_pairwise_Fst_all_scenarios_with_cautions.tsv')

marker_out <- markers
marker_out$interpretation <- ifelse(as.numeric(marker_out$BH_q_four_markers)<0.05,
    'Within-cohort state association; no prospective validity',
    'No BH-adjusted significance within four tested SNPs')
write_tab(marker_out, 'Table_31D_four_marker_associations.tsv')

gates <- data.frame(
  gate=c('Input data readiness','Haploid consensus and site QC',
         'Pairwise PI4K differentiation assessed',
         'Five-state SNP-only geographic assignment',
         'Incremental gain beyond QC supported',
         'Independent geographic validation',
         'Prospective assay validation'),
  status=c('COMPLETE: 282 full; 119 same-study 2020',
           'COMPLETE: 4 strict SNPs in 2020; 18 in full cohort',
           'COMPLETE: sensitivity to marker selection; linked one-gene region',
           'LOW INTERNAL PERFORMANCE: SNP-only BA 26.9%',
           if(lo<=0 && hi>=0) 'UNCERTAIN: conditional interval contains zero' else
             'CONDITIONAL ESTIMATE ONLY; independent assessment needed',
           'NOT DONE','NOT DONE'),
  policy=c('Document archived input checksums and provenance',
           'Do not infer haplotypes from original two-slot coding',
           'Do not interpret as genome-wide Fst',
           'Do not call the panel a validated geographic barcode',
           'Do not turn conditional bootstrap interval into a formal p-value',
           'Use independent collection cohorts with blinded state labels',
           'Evaluate primer specificity, concordance and assay sensitivity')
)
write_tab(gates,'Table_31E_evidence_and_publication_gates.tsv')

state_delta <- 100*(state$mean_combined_recall-state$mean_QC_recall)
report <- c(
 '# PI4K Southwest Nigeria: consolidated analysis evidence (Steps 18–30)',
 '', '## Data and analysis status',
 'Pf8-derived SNP-only PI4K datasets were audited for 282 isolates across five Southwestern Nigerian states.',
 'The principal comparable cohort contains 119 isolates from one contributing study and year (2020).',
 'Four strict SNPs pass the prespecified site-QC filters in 2020. The full cohort has 18 strict SNPs.',
 '', '## Held-out geographic classification',
 sprintf('Repeated nested CV balanced accuracy: SNP-only %.2f%%; QC-only %.2f%%; QC + PI4K %.2f%%.',
         100*snp,100*qc,100*combo),
 sprintf('Mean combined minus QC difference: %.2f percentage points.',100*observed),
 sprintf('Fixed-prediction isolate-bootstrap percentile interval: %.2f to %.2f percentage points.',
         100*lo,100*hi),
 sprintf('Over five correlated CV repeats, PI4K+QC corrected %d QC errors and broke %d correct predictions (net %+d of %d repeated predictions).',
         fixes,breaks,net,119*5),
 'The interval is descriptive conditional on previously fitted models; it is not an independent validation CI or conditional genetic-signal test.',
 '', '## State-specific change in mean recall',
 paste(sprintf('- %s: %+.2f percentage points',state$state,state_delta),collapse='\n'),
 '', '## Differentiation and biological interpretation',
 'Single-gene PI4K differentiation is exploratory and changes with marker policy.',
 'The four marker positions are within one physically linked genomic region.',
 'Some pairs have fewer than two pair-informative SNPs and are not reliably estimable.',
 '', '## Publication claims supported',
 '- Genetically variable PI4K positions and their state frequencies can be reported with qualified evidence.',
 '- Cross-validated classification and its QC comparator can be reported as exploratory internal validation.',
 '- The classifier is not suitable as a validated five-state assignment barcode on current evidence.',
 '', '## Next required independent evidence',
 '- A separately collected cohort from all five states with state-of-infection evidence and travel history.',
 '- Adequate sample numbers per state, balanced assay quality and independent laboratory batches.',
 '- A preregistered locked-marker validation and an expanded multi-locus neutral marker panel.',
 '- No claim of geographic source attribution or clinical deployment without independent validation.',
 '', '## Provenance',
 paste0('Analysis script: R/09_PI4K_final_evidence_consolidation.R. Source directory: ',base),
 'No statistical models were refitted and no original genotypes altered by this summary script.'
)
writeLines(report, file.path(dest,'Step31_consolidated_evidence_report.md'))

if (requireNamespace('ggplot2',quietly=TRUE)) {
  library(ggplot2)
  ggtheme <- theme_classic(base_size=12)+theme(plot.title=element_text(face='bold'),
                                         plot.caption=element_text(size=9,hjust=0))
  order <- c('majority','SNP_nested','SNP_fixed4_posthoc','QC_only','QC_plus_SNP')
  mm <- model_out[match(order,model_out$model), , drop=FALSE]
  mm$model <- factor(mm$model, levels=rev(order))
  p <- ggplot(mm,aes(x=model,y=balanced_accuracy_pct))+
       geom_col(fill='#2766A2',width=.65)+coord_flip()+
       geom_hline(yintercept=20,linetype='dashed',colour='grey45')+
       geom_text(aes(label=sprintf('%.1f%%',balanced_accuracy_pct)),hjust=-.12,size=3.7)+
       scale_y_continuous(limits=c(0,55),expand=c(0,0))+
       labs(title='PI4K geographic classification: nested internal validation',
         subtitle='Five repeated cross-validation partitions of the same 119 isolates',
         x=NULL,y='Mean balanced accuracy (%)',
         caption='Post-hoc four-SNP panel is descriptive only; chance reference is 20%.')+ggtheme
  s <- data.frame(state=state$state,delta_recall_pp=state_delta)
  s$state <- factor(s$state,levels=s$state[order(s$delta_recall_pp)])
  q <- ggplot(s,aes(state,delta_recall_pp))+
       geom_col(fill='#2766A2',width=.65)+coord_flip()+
       geom_hline(yintercept=0,linetype='dashed',colour='grey45')+
       labs(title='Where PI4K changes state-specific recall',
            subtitle='Combined QC+PI4K model minus QC-only, repeated nested CV',
            x=NULL,y='Recall difference (percentage points)',
            caption='Repeated predictions of the same 119 isolates; exploratory internal comparison.')+ggtheme
  save_fig <- function(p, name, width=9,height=5.3) {
    ggplot2::ggsave(file.path(fig,paste0(name,'.pdf')),p,width=width,height=height,
                    units='in',device='pdf',useDingbats=FALSE)
    ggplot2::ggsave(file.path(fig,paste0(name,'.tiff')),p,width=width,height=height,
                    units='in',dpi=600,compression='lzw')
  }
  save_fig(p,'Figure_31A_final_model_comparison')
  save_fig(q,'Figure_31B_recall_change_by_state')
  cat('Created two vector PDFs and two 600-dpi TIFFs.\n')
} else message('ggplot2 unavailable: evidence tables and markdown report still generated')

writeLines(capture.output(sessionInfo()),file.path(dest,'R_sessionInfo.txt'))
cat('STEP 31 COMPLETE: evidence consolidation only; no models refitted.\n')
cat(sprintf('QC+SNP versus QC BA difference: %+.2f pp; conditional interval: %.2f to %.2f pp.\n',
            100*observed,100*lo,100*hi))
cat('FINAL SCIENTIFIC GATE: INDEPENDENT BARCODE VALIDATION NOT COMPLETE.\n')
