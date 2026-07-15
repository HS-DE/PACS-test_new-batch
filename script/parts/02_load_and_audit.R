# -------------------------- Load data --------------------------

study_meta <- standardize_metadata(
  file.path(DATA_DIR, "metadata", "all_sample_metadata.xlsx"), study = TRUE
)
qc_meta <- standardize_metadata(
  file.path(DATA_DIR, "metadata", "all_QC_metadata.xlsx"), study = FALSE
)

study_meta <- study_meta[!study_meta$sample_id %in% qc_meta$sample_id, , drop = FALSE]
study_meta <- study_meta[!study_meta$sample_id %in% OUTLIERS, , drop = FALSE]
qc_meta_use <- qc_meta[qc_meta$sample_role %in% c("internal_qc", "plate_qc"), , drop = FALSE]
all_meta <- rbind(study_meta, qc_meta_use)
if (anyDuplicated(all_meta$sample_id)) stop("Duplicated sample_id after metadata merge.")

protein_path <- file.path(DATA_DIR, "pg_matrix", "pg_matrix.tsv")
protein_raw <- fread(protein_path, data.table = FALSE, check.names = FALSE, na.strings = c("", "NA"))
if (ncol(protein_raw) < 6) stop("Protein matrix does not contain five annotation columns plus samples.")
protein_anno <- protein_raw[, 1:5, drop = FALSE]
protein_mat <- as.matrix(data.frame(lapply(protein_raw[, -(1:5), drop = FALSE], as.numeric), check.names = FALSE))
rownames(protein_mat) <- make.unique(as.character(protein_anno[[1]]))
colnames(protein_mat) <- names(protein_raw)[-(1:5)]
protein_log <- as_log2_protein(protein_mat)

known_ids <- unique(c(study_meta$sample_id, qc_meta$sample_id))
help_log <- read_help_matrix(file.path(DATA_DIR, "HELP", "ALL-HELP.xlsx"), known_ids)

study_ids <- intersect(study_meta$sample_id, intersect(colnames(protein_log), colnames(help_log)))
internal_ids <- intersect(
  qc_meta_use$sample_id[qc_meta_use$sample_role == "internal_qc"],
  intersect(colnames(protein_log), colnames(help_log))
)
plate_qc_ids <- intersect(
  qc_meta_use$sample_id[qc_meta_use$sample_role == "plate_qc"],
  intersect(colnames(protein_log), colnames(help_log))
)

study_meta <- study_meta[match(study_ids, study_meta$sample_id), , drop = FALSE]
internal_meta <- qc_meta_use[match(internal_ids, qc_meta_use$sample_id), , drop = FALSE]
plate_qc_meta <- qc_meta_use[match(plate_qc_ids, qc_meta_use$sample_id), , drop = FALSE]

study_protein <- protein_log[, study_ids, drop = FALSE]
internal_protein <- protein_log[, internal_ids, drop = FALSE]
plate_qc_protein <- protein_log[, plate_qc_ids, drop = FALSE]
study_help <- help_log[, study_ids, drop = FALSE]
internal_help <- help_log[, internal_ids, drop = FALSE]
plate_qc_help <- help_log[, plate_qc_ids, drop = FALSE]

sample_audit <- rbind(
  transform(study_meta, in_protein = sample_id %in% colnames(protein_log), in_help = sample_id %in% colnames(help_log)),
  transform(qc_meta, in_protein = sample_id %in% colnames(protein_log), in_help = sample_id %in% colnames(help_log))
)
sample_audit$protein_detected <- vapply(
  sample_audit$sample_id,
  function(id) if (id %in% colnames(protein_log)) sum(is.finite(protein_log[, id])) else NA_integer_,
  integer(1)
)
sample_audit$help_detected <- vapply(
  sample_audit$sample_id,
  function(id) if (id %in% colnames(help_log)) sum(is.finite(help_log[, id])) else NA_integer_,
  integer(1)
)
fwrite(sample_audit, file.path(OUT_DIR, "audit", "sample_audit.tsv"), sep = "\t")

design_table <- as.data.frame(with(study_meta, table(group2, Plate, Run)))
design_table <- design_table[design_table$Freq > 0, ]
fwrite(design_table, file.path(OUT_DIR, "audit", "group2_plate_run_table.tsv"), sep = "\t")

help_audit <- data.frame(
  HELP = rownames(help_log),
  detect_study = rowMeans(is.finite(study_help)),
  detect_internal = rowMeans(is.finite(internal_help)),
  detect_plate_qc = rowMeans(is.finite(plate_qc_help)),
  mad_study = apply(study_help, 1, mad_or_na),
  mad_internal = apply(internal_help, 1, mad_or_na),
  mad_plate_qc = apply(plate_qc_help, 1, mad_or_na)
)
fwrite(help_audit, file.path(OUT_DIR, "audit", "help_audit.tsv"), sep = "\t")

protein_audit <- data.frame(
  Protein.Group = rownames(study_protein),
  detect_study = rowMeans(is.finite(study_protein)),
  detect_internal = rowMeans(is.finite(internal_protein)),
  detect_plate_qc = rowMeans(is.finite(plate_qc_protein))
)
fwrite(protein_audit, file.path(OUT_DIR, "audit", "protein_audit.tsv"), sep = "\t")

scale_report <- data.frame(
  dataset = c("protein", "HELP"),
  scale = c("log2 transformed from raw positive intensity", attr(help_log, "scale_note"))
)
fwrite(scale_report, file.path(OUT_DIR, "audit", "scale_report.tsv"), sep = "\t")
