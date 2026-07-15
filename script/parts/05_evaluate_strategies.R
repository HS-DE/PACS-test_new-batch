# -------------------------- Evaluation --------------------------

eval_features <- stable_features(raw, study_ids, min_detect = 0.70, max_n = 2500)
raw_group_r2 <- median(partial_r2_term(raw[eval_features, , drop = FALSE], study_meta, "group2"), na.rm = TRUE)

summary_rows <- list()
pca_all <- list()

for (nm in names(strategies)) {
  message("Evaluating ", nm)
  m <- strategies[[nm]]$matrix
  me <- m[eval_features, , drop = FALSE]

  run_r2 <- partial_r2_term(me, study_meta, "Run")
  plate_r2 <- partial_r2_term(me, study_meta, "Plate")
  group_r2 <- partial_r2_term(me, study_meta, "group2")
  pca <- pca_summary(m, study_meta)
  pca_all[[nm]] <- pca

  sample_medians <- apply(m, 2, median_or_na)
  within_group_sd <- tapply(sample_medians, study_meta$group2[match(names(sample_medians), study_meta$sample_id)], sd, na.rm = TRUE)
  shift <- apply(abs(m - raw), 2, median_or_na)

  qc_met <- if (is.null(strategies[[nm]]$qc_matrix)) rep(NA_real_, 3) else control_metrics(strategies[[nm]]$qc_matrix)
  int_met <- if (is.null(strategies[[nm]]$internal_matrix)) rep(NA_real_, 3) else control_metrics(strategies[[nm]]$internal_matrix)

  summary_rows[[nm]] <- data.frame(
    strategy = nm,
    median_run_partial_r2 = median(run_r2, na.rm = TRUE),
    median_plate_partial_r2 = median(plate_r2, na.rm = TRUE),
    median_group_partial_r2 = median(group_r2, na.rm = TRUE),
    group_r2_change_from_raw = median(group_r2, na.rm = TRUE) - raw_group_r2,
    sample_median_sd = sd(sample_medians, na.rm = TRUE),
    mean_within_group_sample_median_sd = mean(within_group_sd, na.rm = TRUE),
    median_absolute_sample_shift = median(shift, na.rm = TRUE),
    p95_absolute_sample_shift = quantile(shift, 0.95, na.rm = TRUE),
    PC1_group_eta2 = pca$eta["PC1_group"],
    PC1_run_eta2 = pca$eta["PC1_run"],
    PC1_plate_eta2 = pca$eta["PC1_plate"],
    QC_median_feature_MAD = qc_met["median_feature_mad"],
    QC_sample_median_SD = qc_met["sample_median_sd"],
    QC_median_pairwise_spearman = qc_met["median_pairwise_spearman"],
    Internal_median_feature_MAD = int_met["median_feature_mad"],
    Internal_sample_median_SD = int_met["sample_median_sd"],
    Internal_median_pairwise_spearman = int_met["median_pairwise_spearman"],
    stringsAsFactors = FALSE
  )

  fwrite(
    data.frame(
      Protein.Group = eval_features,
      Run_partial_R2 = run_r2,
      Plate_partial_R2 = plate_r2,
      Group_partial_R2 = group_r2
    ),
    file.path(OUT_DIR, "metrics", paste0(nm, "_protein_partial_r2.tsv")),
    sep = "\t"
  )

  fwrite(pca$scores, file.path(OUT_DIR, "metrics", paste0(nm, "_pca_scores.tsv")), sep = "\t")

  p <- ggplot(pca$scores, aes(PC1, PC2, colour = group2)) +
    geom_point(size = 2, alpha = 0.75) +
    theme_bw(base_size = 11) +
    labs(
      title = nm,
      subtitle = sprintf("PC1 %.1f%%; PC2 %.1f%%", 100*pca$variance[1], 100*pca$variance[2]),
      colour = "group2"
    )
  ggsave(file.path(OUT_DIR, "plots", paste0(nm, "_PCA.png")), p, width = 7, height = 5, dpi = 160)

  saveRDS(
    list(annotation = protein_anno, matrix = m, metadata = study_meta),
    file.path(OUT_DIR, "matrices", paste0(nm, ".rds")),
    compress = "gzip"
  )
}

strategy_summary <- rbindlist(summary_rows, fill = TRUE)
strategy_summary <- merge(strategy_summary, strategy_info, by = "strategy", all.x = TRUE, sort = FALSE)

rank_low <- function(x) rank(x, na.last = "keep", ties.method = "average")
rank_high <- function(x) rank(-x, na.last = "keep", ties.method = "average")
strategy_summary$screening_score <-
  rank_low(strategy_summary$median_run_partial_r2) +
  rank_low(strategy_summary$median_plate_partial_r2) +
  rank_low(strategy_summary$QC_median_feature_MAD) +
  rank_high(strategy_summary$QC_median_pairwise_spearman) +
  rank_low(abs(strategy_summary$group_r2_change_from_raw)) +
  rank_low(strategy_summary$median_absolute_sample_shift)

strategy_summary <- strategy_summary[order(strategy_summary$screening_score), ]
fwrite(strategy_summary, file.path(OUT_DIR, "metrics", "strategy_summary.tsv"), sep = "\t")
