# Run order and dependencies

These examples document the source workflow in general terms. They will not run until the configured input files exist. Execute commands from the repository root.

**Core analysis stages**

```bash
Rscript R/02_preflight_qc.R
Rscript R/02b_site_selection_audit.R
Rscript R/03_diversity_structure.R
Rscript R/03b_marker_allele_review.R
Rscript R/04b_audited_geographic_differentiation.R
Rscript R/04c_marker_independence_and_FST_robustness.R
```

**Exploratory modelling and follow-up audits**

```bash
Rscript R/05_barcode_nested_cv.R
Rscript R/05a_audit_SNP_vs_QC_model.R
Rscript R/05b_fold_aligned_model_comparison.R
Rscript R/05c_validate_step29_outputs.R
Rscript R/05d_finalize_step29_figures.R
Rscript R/05e_incremental_signal_and_prediction_changes.R
```

**Additional optional reporting steps**

```bash
Rscript R/06_assay_readiness.R
Rscript R/07_finish_audit.R
Rscript R/08_publication_tables.R
Rscript R/09_PI4K_final_evidence_consolidation.R
```

`04_differentiation.R` is retained as original code for transparency; `04b_...` is the more explicit differentiation variant. Consult source preconditions for optional stages. Running one stage may create tabular summaries, graphics, logs and RDS files in the configured research output directory. Those outputs should remain private until disclosure is authorised.

Do not treat successful execution as confirmation that a marker has been externally validated. Train/test leakage checks and predictive model interpretation require an additional methodological review.
