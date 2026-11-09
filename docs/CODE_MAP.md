# Code map

The archived R source files are reproduced without modification. This map describes **general program responsibilities**, not study outcomes. It is not a guarantee that all modules are independently runnable.

| Source file | General responsibility |
| --- | --- |
| `00_config.R` | Input paths, cohort configuration, expected structure and analysis options |
| `01_helpers.R` | File readers, allele coding, metadata loading, statistics and plot export helpers |
| `02_preflight_qc.R` | Input integrity, haploid recoding, preflight sample and site checks |
| `02b_site_selection_audit.R` | Site-filter intersection and allele-count audit |
| `03_diversity_structure.R` | Allele-frequency, diversity, rarefaction and exploratory structure analysis |
| `03b_marker_allele_review.R` | Nucleotide identity, allele-code and within-group frequency review |
| `04_differentiation.R` | Initial differentiation, permutations and associations |
| `04b_audited_geographic_differentiation.R` | More explicit differentiation and marker-policy comparison |
| `04c_marker_independence_and_FST_robustness.R` | Genotype combinations, LD and robustness checks |
| `05_barcode_nested_cv.R` | Exploratory training-fold-selected marker classification |
| `05a_audit_SNP_vs_QC_model.R` | Descriptive classifier-comparator audit |
| `05b_fold_aligned_model_comparison.R` | Matched-fold baseline, QC and genotype model comparisons |
| `05c_validate_step29_outputs.R` | Checks on saved held-out predictions and training metadata |
| `05d_finalize_step29_figures.R` | Generates plots from saved model results without refitting |
| `05e_incremental_signal_and_prediction_changes.R` | Analysis of matched prediction changes and uncertainty |
| `06_assay_readiness.R` | Reference-flank and annotation-overlap preparation |
| `07_finish_audit.R` | Analysis provenance and validation-data templates |
| `08_publication_tables.R` | Format selected model results for reporting |
| `09_PI4K_final_evidence_consolidation.R` | Consolidate saved analyses into technical summaries |

Some modules share state through intermediate serialized R objects and result tables. Downstream steps require successful execution of their dependencies. This is a research workflow, not a supported stand-alone software package.
