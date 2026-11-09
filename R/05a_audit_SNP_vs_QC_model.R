#!/usr/bin/env Rscript
# Step 28: Post-hoc audit of saved Step 27 outputs (not refitting classifiers).
# Reproducible, dependency-light: base R + ggplot2 for figures.
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1]], mustWork = TRUE) else "/mnt/c/PI4K_Nigeria_Study"
base <- file.path(root, "results", "Southwest5", "R_Analysis_v1")
model <- file.path(base, "model")
tables <- file.path(base, "tables")
figdir <- file.path(base, "figures")
dir.create(tables, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
readtsv <- function(name) {
  path <- file.path(model, name)
  if (!file.exists(path)) stop("Missing input: ", path)
  read.delim(path, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE, quote = '"')
}
need <- function(df, cols, name) {
  missing <- setdiff(cols, names(df))
  if (length(missing)) stop(name, ": missing column(s): ", paste(missing, collapse = ", "))
}
writeout <- function(df, name) {
  path <- file.path(tables, name)
  write.table(df, path, sep = "\t", quote = TRUE, row.names = FALSE, na = "NA")
  message("Saved: ", path)
}
met <- readtsv("Table_5_nested_CV_metrics_by_repeat.tsv")
need(met, c("repeat_id", "balanced_accuracy", "macro_F1", "overall_accuracy", "coverage",
            "QC_only_balanced_accuracy", "QC_only_macro_F1", "QC_only_overall_accuracy", "QC_only_coverage"), "CV metrics")
for (nm in setdiff(names(met), "repeat_id")) met[[nm]] <- as.numeric(met[[nm]])
if (any(!is.finite(met$balanced_accuracy)) || any(!is.finite(met$QC_only_balanced_accuracy))) stop("Non-finite balanced accuracy")
if (any(duplicated(met$repeat_id))) stop("Duplicate repeat_id")
met$delta_balanced_accuracy_SNP_minus_QC <- met$balanced_accuracy - met$QC_only_balanced_accuracy
met$delta_macro_F1_SNP_minus_QC <- met$macro_F1 - met$QC_only_macro_F1
met$delta_coverage_SNP_minus_QC <- met$coverage - met$QC_only_coverage
writeout(met, "Table_S36_repeated_CV_SNP_vs_QC_paired_descriptive.tsv")

cm <- readtsv("Table_S17_first_repeat_confusion_matrix.tsv")
rec <- readtsv("Table_S18_first_repeat_per_state_recall.tsv")
need(cm, "observed_state", "confusion matrix")
need(rec, c("state", "n", "correct", "recall"), "recall table")
states <- as.character(cm$observed_state)
if (anyDuplicated(states)) stop("Duplicate state in confusion table")
if (!setequal(states, names(cm)[-1])) stop("Confusion matrix column states differ from row states")
if (!setequal(states, rec$state)) stop("Recall states do not match confusion matrix")
rec <- rec[match(states, rec$state), , drop = FALSE]
mat <- as.matrix(cm[, states, drop = FALSE])
storage.mode(mat) <- "numeric"
if (any(!is.finite(mat)) || any(mat < 0)) stop("Invalid confusion counts")
row.names(mat) <- states
n_obs <- as.numeric(rec$n)
n_pred <- rowSums(mat)
correct <- diag(mat)
if (any(n_obs < n_pred)) stop("Predictions exceed observed state total")
if (any(correct != as.numeric(rec$correct))) stop("Recall and confusion diagonal do not agree")
pred_col <- colSums(mat)
comparison <- data.frame(
  state = states,
  n_true = n_obs,
  n_predictions_made = n_pred,
  n_abstentions = n_obs - n_pred,
  n_correct = correct,
  recall_with_abstention_as_failure = correct / n_obs,
  recall_conditional_on_prediction = ifelse(n_pred > 0, correct / n_pred, NA_real_),
  n_predicted_as_this_state = pred_col[states],
  precision_when_predicted = ifelse(pred_col[states] > 0, correct / pred_col[states], NA_real_),
  stringsAsFactors = FALSE
)
writeout(comparison, "Table_S37_first_repeat_prediction_abstention_audit.tsv")

summ <- readtsv("Table_5_barcode_internal_validation_summary.tsv")
need(summ, c("model", "whole_pipeline_permutation_p", "n_whole_pipeline_null_permutations"), "main model summary")
summary <- data.frame(
  metric = c("n_repeats", "SNP_balanced_accuracy_mean", "QC_balanced_accuracy_mean", "SNP_minus_QC_balanced_accuracy_mean",
             "SNP_macro_F1_mean", "QC_macro_F1_mean", "SNP_minus_QC_macro_F1_mean",
             "SNP_prediction_coverage_mean", "QC_prediction_coverage_mean", "SNP_minus_QC_coverage_mean",
             "SNP_balanced_accuracy_better_than_QC_n_repeats", "first_repeat_prediction_n", "first_repeat_total_n",
             "first_repeat_abstentions_n", "first_repeat_states_never_predicted_n",
             "published_single_repeat_pipeline_permutation_p", "pipeline_null_permutations_n"),
  value = as.character(c(nrow(met), mean(met$balanced_accuracy), mean(met$QC_only_balanced_accuracy),
                         mean(met$delta_balanced_accuracy_SNP_minus_QC), mean(met$macro_F1), mean(met$QC_only_macro_F1),
                         mean(met$delta_macro_F1_SNP_minus_QC), mean(met$coverage), mean(met$QC_only_coverage),
                         mean(met$delta_coverage_SNP_minus_QC), sum(met$delta_balanced_accuracy_SNP_minus_QC > 0),
                         sum(n_pred), sum(n_obs), sum(n_obs - n_pred), sum(pred_col == 0),
                         summ$whole_pipeline_permutation_p[1], summ$n_whole_pipeline_null_permutations[1])),
  note = "Descriptive audit only; repeated CV runs are correlated, no paired inferential p-value",
  stringsAsFactors = FALSE
)
writeout(summary, "Table_S38_step27_model_comparison_audit_summary.tsv")

if (requireNamespace("ggplot2", quietly = TRUE)) {
  gg <- asNamespace("ggplot2")
  # Use fixed model labels; no validation data are modified.
  plotdf <- rbind(
    data.frame(repeat_id = as.character(met$repeat_id), model = "PI4K SNPs", balanced_accuracy = met$balanced_accuracy),
    data.frame(repeat_id = as.character(met$repeat_id), model = "QC-only", balanced_accuracy = met$QC_only_balanced_accuracy)
  )
  plotdf$model <- factor(plotdf$model, levels = c("PI4K SNPs", "QC-only"))
  plotdf$repeat_id <- factor(plotdf$repeat_id, levels = as.character(met$repeat_id))
  p1 <- ggplot2::ggplot(plotdf, ggplot2::aes(x = repeat_id, y = balanced_accuracy, color = model, group = model)) +
    ggplot2::geom_line(linewidth = 0.75) + ggplot2::geom_point(size = 2.8) +
    ggplot2::geom_hline(yintercept = 0.2, linetype = "dashed", color = "grey60") +
    ggplot2::scale_y_continuous(labels = function(x) sprintf("%.0f%%", 100 * x), limits = c(0, max(0.5, plotdf$balanced_accuracy) * 1.1)) +
    ggplot2::labs(title = "PI4K SNPs versus quality-only features", subtitle = "Five nested-CV repeats; descriptive comparison, not independent replication",
                    x = "Outer CV repeat", y = "Balanced accuracy", color = NULL,
                    caption = "Dashed line: five-state chance reference (20%). Models may differ in prediction coverage.") +
    ggplot2::theme_minimal(base_size = 12) + ggplot2::theme(legend.position = "top", panel.grid.minor = ggplot2::element_blank())
  p2df <- comparison
  p2df$state <- factor(p2df$state, levels = states)
  p2 <- ggplot2::ggplot(p2df, ggplot2::aes(x = state, y = recall_with_abstention_as_failure)) +
    ggplot2::geom_col(fill = "#356B84", width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = paste0(n_correct, "/", n_true)), vjust = -0.5, size = 3.7) +
    ggplot2::scale_y_continuous(labels = function(x) sprintf("%.0f%%", 100 * x), limits = c(0, max(0.7, p2df$recall_with_abstention_as_failure * 1.12))) +
    ggplot2::labs(title = "State-wise recall, first outer-CV repeat", subtitle = "Missing predictions counted as incorrect",
                    x = "State of sample collection", y = "Recall",
                    caption = "Only first-repeat results; not a population-level performance estimate.") +
    ggplot2::theme_minimal(base_size = 12) + ggplot2::theme(panel.grid.major.x = ggplot2::element_blank())
  savefig <- function(p, stem, w, h) {
    pdf <- file.path(figdir, paste0(stem, ".pdf"))
    tif <- file.path(figdir, paste0(stem, ".tiff"))
    ggplot2::ggsave(pdf, p, width = w, height = h, units = "in", device = grDevices::cairo_pdf)
    ggplot2::ggsave(tif, p, width = w, height = h, units = "in", device = "tiff", dpi = 600, compression = "lzw")
    message("Saved: ", pdf, " and ", tif)
  }
  savefig(p1, "Figure_15_SNP_vs_QC_repeated_balanced_accuracy_AUDIT", 8, 4.8)
  savefig(p2, "Figure_16_first_repeat_state_recall_AUDIT", 8, 4.8)
} else {
  message("ggplot2 unavailable: tables exported; figures not generated")
}
message(sprintf("Step 28 audit completed: mean SNP BA %.4f; QC-only BA %.4f; delta %.4f", mean(met$balanced_accuracy), mean(met$QC_only_balanced_accuracy), mean(met$delta_balanced_accuracy_SNP_minus_QC)))
message("CAUTION: Permutation p-value tests label independence for a specified score; it does not test superiority to QC-only or clinical utility.")
