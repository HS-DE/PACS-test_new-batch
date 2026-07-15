# -------------------------- HELP sets --------------------------

all_analysis_ids <- c(study_ids, internal_ids, plate_qc_ids)
complete_help <- rownames(help_log)[rowSums(is.finite(help_log[, all_analysis_ids, drop = FALSE])) == length(all_analysis_ids)]

eligible_help <- help_audit$HELP[
  help_audit$detect_study >= 0.80 &
  help_audit$detect_internal >= 0.67 &
  help_audit$detect_plate_qc >= 0.67
]
if (length(eligible_help) < 10) {
  eligible_help <- help_audit$HELP[
    help_audit$detect_study >= 0.60 &
    help_audit$detect_internal >= 0.50 &
    help_audit$detect_plate_qc >= 0.50
  ]
}
if (length(complete_help) < 5) complete_help <- eligible_help

help_weights <- 1 / (help_audit$mad_study^2 + 0.05)
names(help_weights) <- help_audit$HELP
finite_w <- help_weights[is.finite(help_weights)]
if (length(finite_w)) {
  cap <- quantile(finite_w, c(0.05, 0.95), na.rm = TRUE)
  help_weights <- pmin(pmax(help_weights, cap[1]), cap[2])
}

help_sets <- data.frame(
  HELP = rownames(help_log),
  complete = rownames(help_log) %in% complete_help,
  eligible_available = rownames(help_log) %in% eligible_help
)
fwrite(help_sets, file.path(OUT_DIR, "audit", "help_sets.tsv"), sep = "\t")

# -------------------------- Factors --------------------------

run_internal_complete_sample <- feature_deviation_factors(
  internal_help, internal_ids, complete_help, min_features = 5
)
run_complete <- aggregate_factor(run_internal_complete_sample, internal_meta, "Run")

run_internal_available_sample <- feature_deviation_factors(
  internal_help, internal_ids, eligible_help, min_features = 5, weights = help_weights
)
run_available <- aggregate_factor(run_internal_available_sample, internal_meta, "Run")

internal_protein_features <- stable_features(internal_protein, internal_ids, 0.75, 2000)
run_internal_protein_sample <- feature_deviation_factors(
  internal_protein, internal_ids, internal_protein_features, min_features = 50
)
run_protein <- aggregate_factor(run_internal_protein_sample, internal_meta, "Run")
dual_run <- merge(run_available, run_protein, by = "batch", all = TRUE, suffixes = c("_help", "_protein"))
dual_run$k <- rowMeans(dual_run[, c("k_help", "k_protein")], na.rm = TRUE)
dual_run$k <- dual_run$k - median(dual_run$k, na.rm = TRUE)
dual_run <- dual_run[, c("batch", "k")]

k_run_complete_study <- factor_for_samples(study_meta, run_complete, "Run")
k_run_complete_qc <- factor_for_samples(plate_qc_meta, run_complete, "Run")
k_run_complete_internal <- factor_for_samples(internal_meta, run_complete, "Run")
k_run_available_study <- factor_for_samples(study_meta, run_available, "Run")
k_run_available_qc <- factor_for_samples(plate_qc_meta, run_available, "Run")
k_run_available_internal <- factor_for_samples(internal_meta, run_available, "Run")
k_dual_study <- factor_for_samples(study_meta, dual_run, "Run")
k_dual_qc <- factor_for_samples(plate_qc_meta, dual_run, "Run")
k_dual_internal <- factor_for_samples(internal_meta, dual_run, "Run")

qc_run_complete <- apply_named_factor(plate_qc_protein, k_run_complete_qc)
qc_features <- stable_features(qc_run_complete, plate_qc_ids, 0.67, 3000)
plate_qc_sample_complete <- feature_deviation_factors(
  qc_run_complete, plate_qc_ids, qc_features, min_features = 50
)
plate_complete <- aggregate_factor(plate_qc_sample_complete, plate_qc_meta, "Plate")

qc_run_available <- apply_named_factor(plate_qc_protein, k_run_available_qc)
qc_features_available <- stable_features(qc_run_available, plate_qc_ids, 0.67, 3000)
plate_qc_sample_available <- feature_deviation_factors(
  qc_run_available, plate_qc_ids, qc_features_available, min_features = 50
)
plate_available <- aggregate_factor(plate_qc_sample_available, plate_qc_meta, "Plate")

qc_run_dual <- apply_named_factor(plate_qc_protein, k_dual_qc)
qc_features_dual <- stable_features(qc_run_dual, plate_qc_ids, 0.67, 3000)
plate_qc_sample_dual <- feature_deviation_factors(
  qc_run_dual, plate_qc_ids, qc_features_dual, min_features = 50
)
plate_dual <- aggregate_factor(plate_qc_sample_dual, plate_qc_meta, "Plate")

k_plate_study <- factor_for_samples(study_meta, plate_complete, "Plate")
k_plate_qc <- factor_for_samples(plate_qc_meta, plate_complete, "Plate")
k_plate_available_study <- factor_for_samples(study_meta, plate_available, "Plate")
k_plate_available_qc <- factor_for_samples(plate_qc_meta, plate_available, "Plate")
k_plate_dual_study <- factor_for_samples(study_meta, plate_dual, "Plate")
k_plate_dual_qc <- factor_for_samples(plate_qc_meta, plate_dual, "Plate")

study_help_run_complete <- apply_named_factor(study_help, k_run_complete_study)
study_help_run_available <- apply_named_factor(study_help, k_run_available_study)
study_help_run_dual <- apply_named_factor(study_help, k_dual_study)

group2_named <- study_meta$group2
names(group2_named) <- study_meta$sample_id

sample_complete_global <- feature_deviation_factors(
  study_help_run_complete, study_ids, complete_help, min_features = 5
)
sample_complete_group <- feature_deviation_factors(
  study_help_run_complete, study_ids, complete_help,
  groups = group2_named, min_features = 5
)
sample_available_global <- feature_deviation_factors(
  study_help_run_available, study_ids, eligible_help,
  min_features = 8, weights = help_weights
)
sample_dual_available <- feature_deviation_factors(
  study_help_run_dual, study_ids, eligible_help,
  min_features = 8, weights = help_weights
)

k_sample_complete_global <- setNames(sample_complete_global$k, sample_complete_global$sample_id)
k_sample_complete_group <- setNames(sample_complete_group$k, sample_complete_group$sample_id)
k_sample_available_global <- setNames(sample_available_global$k, sample_available_global$sample_id)
k_sample_dual_available <- setNames(sample_dual_available$k, sample_dual_available$sample_id)

legacy_internal_sample <- column_median_factors(internal_help, internal_ids, min_features = 5)
legacy_run <- aggregate_factor(legacy_internal_sample, internal_meta, "Run")
legacy_qc_help_sample <- column_median_factors(plate_qc_help, plate_qc_ids, min_features = 5)
legacy_qc_help_plate <- aggregate_factor(legacy_qc_help_sample, plate_qc_meta, "Plate")
legacy_qc1 <- apply_named_factor(plate_qc_protein, factor_for_samples(plate_qc_meta, legacy_run, "Run"))
legacy_qc21 <- apply_named_factor(legacy_qc1, factor_for_samples(plate_qc_meta, legacy_qc_help_plate, "Plate"))
legacy_qc_protein_sample <- column_median_factors(legacy_qc21, plate_qc_ids, min_features = 50)
legacy_plate <- aggregate_factor(legacy_qc_protein_sample, plate_qc_meta, "Plate")
legacy_study_help_sample <- column_median_factors(
  study_help, study_ids, groups = group2_named, min_features = 5
)

factor_export <- rbind(
  transform(run_complete, factor_type = "run_complete"),
  transform(run_available, factor_type = "run_available"),
  transform(dual_run, factor_type = "run_dual"),
  transform(plate_complete, factor_type = "plate_complete"),
  transform(plate_available, factor_type = "plate_available"),
  transform(plate_dual, factor_type = "plate_dual"),
  transform(legacy_run, factor_type = "legacy_run"),
  transform(legacy_plate, factor_type = "legacy_plate")
)
fwrite(factor_export, file.path(OUT_DIR, "factors", "batch_factors.tsv"), sep = "\t")

sample_factor_export <- rbind(
  transform(sample_complete_global, factor_type = "sample_complete_global"),
  transform(sample_complete_group, factor_type = "sample_complete_group"),
  transform(sample_available_global, factor_type = "sample_available_global"),
  transform(sample_dual_available, factor_type = "sample_dual_available"),
  transform(legacy_study_help_sample[, c("sample_id", "k", "n_features")], factor_type = "legacy_group_help")
)
fwrite(sample_factor_export, file.path(OUT_DIR, "factors", "sample_factors.tsv"), sep = "\t")
