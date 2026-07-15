# -------------------------- Strategies --------------------------

strategies <- list()
strategy_info <- data.frame(
  strategy = character(), family = character(), uses_group_in_preprocessing = logical(),
  description = character(), stringsAsFactors = FALSE
)
add_strategy <- function(name, matrix, family, uses_group, description, qc_matrix = NULL, internal_matrix = NULL) {
  strategies[[name]] <<- list(matrix = matrix, qc_matrix = qc_matrix, internal_matrix = internal_matrix)
  strategy_info <<- rbind(strategy_info, data.frame(
    strategy = name, family = family,
    uses_group_in_preprocessing = uses_group,
    description = description, stringsAsFactors = FALSE
  ))
}

raw <- study_protein
run_only <- apply_named_factor(raw, k_run_complete_study)
run_plate <- apply_named_factor(run_only, k_plate_study)
seq_complete <- apply_named_factor(run_plate, k_sample_complete_global)
seq_complete_gm <- group_median_normalize(seq_complete, study_meta, "group2")
seq_group_help <- apply_named_factor(run_plate, k_sample_complete_group)
seq_group_help_gm <- group_median_normalize(seq_group_help, study_meta, "group2")
seq_available <- apply_named_factor(
  apply_named_factor(apply_named_factor(raw, k_run_available_study), k_plate_available_study),
  k_sample_available_global
)
dual_base <- apply_named_factor(apply_named_factor(raw, k_dual_study), k_plate_dual_study)
dual_available <- apply_named_factor(dual_base, k_sample_dual_available)

legacy <- apply_named_factor(raw, factor_for_samples(study_meta, legacy_run, "Run"))
legacy <- apply_named_factor(legacy, factor_for_samples(study_meta, legacy_plate, "Plate"))
legacy_k_sample <- setNames(legacy_study_help_sample$k, legacy_study_help_sample$sample_id)
legacy <- apply_named_factor(legacy, legacy_k_sample)
legacy_gm <- group_median_normalize(legacy, study_meta, "group2")$matrix

help_pc <- make_help_pcs(study_help_run_available, study_ids, eligible_help, n_pc = 2)
if (!is.null(help_pc)) {
  fwrite(help_pc$scores, file.path(OUT_DIR, "factors", "help_pc_scores.tsv"), sep = "\t")
  fwrite(data.frame(PC = seq_along(help_pc$variance), variance = help_pc$variance),
         file.path(OUT_DIR, "factors", "help_pc_variance.tsv"), sep = "\t")
}

batch_model <- protected_adjust(raw, study_meta)
batch_help_model <- protected_adjust(raw, study_meta, if (is.null(help_pc)) NULL else help_pc$scores)
conditional_model <- protected_adjust(
  raw, study_meta, if (is.null(help_pc)) NULL else help_pc$scores,
  conditional_threshold = 0.05
)

qc_raw <- plate_qc_protein
internal_raw <- internal_protein
qc_run <- apply_named_factor(qc_raw, k_run_complete_qc)
internal_run <- apply_named_factor(internal_raw, k_run_complete_internal)
qc_run_plate <- apply_named_factor(qc_run, k_plate_qc)
qc_available <- apply_named_factor(apply_named_factor(qc_raw, k_run_available_qc), k_plate_available_qc)
qc_dual <- apply_named_factor(apply_named_factor(qc_raw, k_dual_qc), k_plate_dual_qc)
internal_dual <- apply_named_factor(internal_raw, k_dual_internal)

add_strategy("S00_raw", raw, "baseline", FALSE, "Raw log2 protein data", qc_raw, internal_raw)
add_strategy("S01_legacy_group2", legacy_gm, "legacy", TRUE, "Legacy-like independent Run/Plate plus group2 HELP and group2 median", legacy_qc21, internal_run)
add_strategy("S02_run_only", run_only, "scalar", FALSE, "Internal HELP Run correction only", qc_run, internal_run)
add_strategy("S03_minimal_run_plate", run_plate, "scalar", FALSE, "Run plus QC endogenous-protein Plate correction", qc_run_plate, internal_run)
add_strategy("S04_sequential_complete", seq_complete, "scalar", FALSE, "Run + Plate + complete-HELP sample correction", qc_run_plate, internal_run)
add_strategy("S05_sequential_complete_group_median", seq_complete_gm$matrix, "scalar", TRUE, "S04 plus group2 median normalization", qc_run_plate, internal_run)
add_strategy("S06_group_help", seq_group_help, "scalar", TRUE, "Run + Plate + group2-specific HELP correction", qc_run_plate, internal_run)
add_strategy("S07_group_help_group_median", seq_group_help_gm$matrix, "scalar", TRUE, "S06 plus group2 median normalization", qc_run_plate, internal_run)
add_strategy("S08_available_help", seq_available, "scalar", FALSE, "Missing-tolerant weighted HELP correction", qc_available, apply_named_factor(internal_raw, k_run_available_internal))
add_strategy("S09_dual_anchor", dual_available, "scalar", FALSE, "Dual Internal HELP/protein Run anchor + Plate + available HELP", qc_dual, internal_dual)
add_strategy("S10_batch_model", batch_model$matrix, "model", FALSE, "Protein-specific protected group2 + Run + Plate model")
add_strategy("S11_batch_help_pc_model", batch_help_model$matrix, "model", FALSE, "Protein-specific protected model with HELP PCs")
add_strategy("S12_conditional_model", conditional_model$matrix, "model", FALSE, "Apply model correction only when technical partial R2 >= 0.05")

fwrite(strategy_info, file.path(OUT_DIR, "strategy_manifest.tsv"), sep = "\t")
