suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)
set.seed(20260714)

ROOT <- normalizePath(".", mustWork = TRUE)
DATA_DIR <- file.path(ROOT, "data")
OUT_DIR <- file.path(ROOT, "Results")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
for (d in c("audit", "factors", "metrics", "plots", "matrices")) {
  dir.create(file.path(OUT_DIR, d), recursive = TRUE, showWarnings = FALSE)
}

OUTLIERS <- c(
  "WJ-1-5","WJ-1-7","WJ-1-19","WJ-1-21","WJ-1-44","WJ-1-67",
  "WJ-1-90","WJ-1-101","WJ-1-102","WJ-1-103","WJ-1-104",
  "WJ-1-105","WJ-1-112","WJ-1-146","WJ-1-148","WJ-1-157",
  "WJ-1-171","WJ-1-176","WJ-1-177","WJ-1-178","WJ-1-187",
  "WJ-1-203","WJ-1-208","WJ-1-209","WJ-1-212","WJ-1-213",
  "WJ-1-43","WJ-1-59","WJ-1-83"
)

clean_name <- function(x) {
  x <- tolower(trimws(as.character(x)))
  gsub("[^a-z0-9]+", "_", x)
}

pick_col <- function(df, candidates, required = TRUE) {
  nm <- clean_name(names(df))
  cand <- clean_name(candidates)
  hit <- which(nm %in% cand)
  if (length(hit)) return(names(df)[hit[1]])
  for (z in cand) {
    hit <- grep(z, nm, fixed = TRUE)
    if (length(hit)) return(names(df)[hit[1]])
  }
  if (required) stop("Cannot find column: ", paste(candidates, collapse = ", "))
  NULL
}

safe_char <- function(df, col, default = NA_character_) {
  if (is.null(col)) return(rep(default, nrow(df)))
  as.character(df[[col]])
}

standardize_metadata <- function(path, study = FALSE) {
  x <- as.data.frame(read_excel(path, sheet = 1, .name_repair = "minimal"))
  sid <- pick_col(x, c("sample_id", "sampleid", "sample id"))
  sname <- pick_col(x, c("sample_name", "samplename", "sample name"), FALSE)
  group1 <- pick_col(x, c("group1", "group_1", "group"), FALSE)
  group2 <- pick_col(x, c("group2", "group_2"), FALSE)
  plate <- pick_col(x, c("plate", "plate_id", "batch"), FALSE)
  run <- pick_col(x, c("run", "run_id", "instrument_run"), FALSE)

  out <- data.frame(
    sample_id = trimws(safe_char(x, sid)),
    sample_name = trimws(safe_char(x, sname)),
    group1 = trimws(safe_char(x, group1)),
    group2 = trimws(safe_char(x, group2)),
    Plate = trimws(safe_char(x, plate)),
    Run = trimws(safe_char(x, run)),
    stringsAsFactors = FALSE
  )

  if (study) {
    out$sample_role <- "study"
  } else {
    all_text <- apply(x, 1, function(z) paste(tolower(as.character(z)), collapse = " | "))
    out$sample_role <- ifelse(
      grepl("blank", all_text), "blank",
      ifelse(
        grepl("neat|qc2", all_text), "neat",
        ifelse(
          grepl("internal|calibrator|qc3", all_text), "internal_qc",
          ifelse(grepl("qc1|plate.?qc", all_text), "plate_qc", "other_control")
        )
      )
    )
  }

  out$group2[is.na(out$group2) | !nzchar(out$group2)] <-
    out$group1[is.na(out$group2) | !nzchar(out$group2)]
  out
}

read_help_matrix <- function(path, known_ids) {
  x <- as.data.frame(read_excel(path, sheet = 1, .name_repair = "minimal"))
  if (ncol(x) < 2) stop("HELP file has fewer than two columns.")

  col_overlap <- sum(names(x) %in% known_ids)
  row_overlap <- sum(as.character(x[[1]]) %in% known_ids)

  if (col_overlap >= row_overlap) {
    feature_id <- as.character(x[[1]])
    m <- as.matrix(data.frame(lapply(x[-1], as.numeric), check.names = FALSE))
    rownames(m) <- make.unique(feature_id)
    colnames(m) <- names(x)[-1]
  } else {
    sample_id <- as.character(x[[1]])
    m <- t(as.matrix(data.frame(lapply(x[-1], as.numeric), check.names = FALSE)))
    rownames(m) <- make.unique(names(x)[-1])
    colnames(m) <- sample_id
  }

  m[!is.finite(m) | m <= 0] <- NA_real_
  med <- median(as.numeric(m[is.finite(m)]), na.rm = TRUE)
  if (!is.finite(med)) stop("No finite HELP values found.")

  if (med > 100) {
    m <- log2(m)
    attr(m, "scale_note") <- paste0("log2 transformed automatically; raw median=", signif(med, 5))
  } else {
    attr(m, "scale_note") <- paste0("treated as already log2; median=", signif(med, 5))
  }
  m
}

as_log2_protein <- function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  x[!is.finite(x) | x <= 0] <- NA_real_
  log2(x)
}

median_or_na <- function(x) {
  if (!any(is.finite(x))) return(NA_real_)
  median(x, na.rm = TRUE)
}

mad_or_na <- function(x) {
  if (sum(is.finite(x)) < 2) return(NA_real_)
  mad(x, na.rm = TRUE, constant = 1)
}

weighted_median <- function(x, w = NULL) {
  ok <- is.finite(x)
  x <- x[ok]
  if (!length(x)) return(NA_real_)
  if (is.null(w)) return(median(x))
  w <- w[ok]
  ok2 <- is.finite(w) & w > 0
  x <- x[ok2]; w <- w[ok2]
  if (!length(x)) return(NA_real_)
  o <- order(x)
  x <- x[o]; w <- w[o] / sum(w)
  x[which(cumsum(w) >= 0.5)[1]]
}

feature_deviation_factors <- function(mat, ids, features = rownames(mat), groups = NULL, min_features = 5L, weights = NULL) {
  ids <- intersect(ids, colnames(mat))
  features <- intersect(features, rownames(mat))
  sub <- mat[features, ids, drop = FALSE]
  result <- data.frame(sample_id = ids, k = NA_real_, n_features = 0L)

  if (is.null(groups)) {
    group_vector <- rep("__all__", length(ids)); names(group_vector) <- ids
  } else {
    group_vector <- as.character(groups[ids])
    group_vector[is.na(group_vector) | !nzchar(group_vector)] <- "__missing__"
  }

  for (g in unique(group_vector)) {
    gid <- names(group_vector)[group_vector == g]
    ref <- apply(sub[, gid, drop = FALSE], 1, median_or_na)
    for (id in gid) {
      d <- sub[, id] - ref
      ok <- is.finite(d)
      result$n_features[result$sample_id == id] <- sum(ok)
      if (sum(ok) >= min_features) {
        result$k[result$sample_id == id] <- -weighted_median(d, if (is.null(weights)) NULL else weights[features])
      }
    }
  }
  result
}

column_median_factors <- function(mat, ids, groups = NULL, min_features = 5L) {
  ids <- intersect(ids, colnames(mat))
  med <- apply(mat[, ids, drop = FALSE], 2, median_or_na)
  out <- data.frame(sample_id = ids, median = med, k = NA_real_)
  if (is.null(groups)) {
    out$k <- median(med, na.rm = TRUE) - out$median
  } else {
    gv <- as.character(groups[ids])
    for (g in unique(gv)) {
      idx <- which(gv == g)
      out$k[idx] <- median(med[idx], na.rm = TRUE) - med[idx]
    }
  }
  out$n_features <- colSums(is.finite(mat[, ids, drop = FALSE]))
  out$k[out$n_features < min_features] <- NA_real_
  out
}

aggregate_factor <- function(sample_factor, meta, group_col, value_col = "k") {
  z <- merge(sample_factor, meta[, c("sample_id", group_col), drop = FALSE], by = "sample_id", all.x = TRUE, sort = FALSE)
  names(z)[names(z) == group_col] <- "batch"
  z <- z[!is.na(z$batch) & nzchar(z$batch), , drop = FALSE]
  ans <- aggregate(z[[value_col]], list(batch = z$batch), median, na.rm = TRUE)
  names(ans)[2] <- "k"
  ans$k <- ans$k - median(ans$k, na.rm = TRUE)
  ans
}

factor_for_samples <- function(meta, batch_factor, batch_col) {
  k <- batch_factor$k[match(as.character(meta[[batch_col]]), batch_factor$batch)]
  names(k) <- meta$sample_id
  k
}

apply_named_factor <- function(mat, k) {
  z <- k[colnames(mat)]
  z[!is.finite(z)] <- 0
  sweep(mat, 2, z, "+")
}

group_median_normalize <- function(mat, meta, group_col = "group2") {
  med <- apply(mat, 2, median_or_na)
  g <- as.character(meta[[group_col]]); names(g) <- meta$sample_id
  k <- rep(NA_real_, ncol(mat)); names(k) <- colnames(mat)
  for (lev in unique(g[colnames(mat)])) {
    ids <- colnames(mat)[g[colnames(mat)] == lev]
    k[ids] <- median(med[ids], na.rm = TRUE) - med[ids]
  }
  list(matrix = apply_named_factor(mat, k), factor = k)
}

stable_features <- function(mat, ids, min_detect = 0.7, max_n = Inf) {
  ids <- intersect(ids, colnames(mat))
  det <- rowMeans(is.finite(mat[, ids, drop = FALSE]))
  variability <- apply(mat[, ids, drop = FALSE], 1, mad_or_na)
  keep <- which(det >= min_detect & is.finite(variability))
  if (length(keep) > max_n) keep <- keep[order(variability[keep])][seq_len(max_n)]
  rownames(mat)[keep]
}

eta2_factor <- function(x, g) {
  ok <- is.finite(x) & !is.na(g)
  x <- x[ok]; g <- as.factor(g[ok])
  if (length(x) < 3 || nlevels(g) < 2) return(NA_real_)
  grand <- mean(x); ss_total <- sum((x - grand)^2)
  if (ss_total <= 0) return(0)
  means <- tapply(x, g, mean); ns <- table(g)
  sum(ns * (means - grand)^2) / ss_total
}

make_help_pcs <- function(help_mat, ids, features, n_pc = 2L) {
  ids <- intersect(ids, colnames(help_mat)); features <- intersect(features, rownames(help_mat))
  x <- t(help_mat[features, ids, drop = FALSE])
  for (j in seq_len(ncol(x))) x[!is.finite(x[, j]), j] <- median(x[, j], na.rm = TRUE)
  x <- x[, apply(x, 2, sd, na.rm = TRUE) > 0, drop = FALSE]
  if (ncol(x) < 2) return(NULL)
  pc <- prcomp(x, center = TRUE, scale. = TRUE)
  n_pc <- min(n_pc, ncol(pc$x))
  scores <- as.data.frame(pc$x[, seq_len(n_pc), drop = FALSE]); names(scores) <- paste0("HELP_PC", seq_len(n_pc)); scores$sample_id <- rownames(scores)
  list(scores = scores, variance = pc$sdev^2 / sum(pc$sdev^2), loadings = pc$rotation[, seq_len(n_pc), drop = FALSE])
}

protected_adjust <- function(mat, meta, continuous = NULL, conditional_threshold = NULL) {
  ids <- intersect(colnames(mat), meta$sample_id)
  meta2 <- meta[match(ids, meta$sample_id), , drop = FALSE]
  ymat <- mat[, ids, drop = FALSE]
  design_df <- data.frame(group2 = factor(meta2$group2), Run = factor(meta2$Run), Plate = factor(meta2$Plate))
  if (!is.null(continuous)) {
    continuous <- continuous[match(ids, continuous$sample_id), , drop = FALSE]
    for (nm in setdiff(names(continuous), "sample_id")) design_df[[nm]] <- as.numeric(continuous[[nm]])
  }
  rhs <- c("group2", "Run", "Plate", setdiff(names(design_df), c("group2", "Run", "Plate")))
  Xfull_all <- model.matrix(as.formula(paste("~", paste(rhs, collapse = " + "))), design_df)
  Xprot_all <- model.matrix(~ group2, design_df)
  tech_cols <- setdiff(colnames(Xfull_all), colnames(Xprot_all))
  corrected <- ymat
  partial_r2 <- rep(NA_real_, nrow(ymat)); names(partial_r2) <- rownames(ymat)
  applied <- rep(FALSE, nrow(ymat))
  for (j in seq_len(nrow(ymat))) {
    y <- ymat[j, ]; ok <- is.finite(y) & apply(is.finite(Xfull_all), 1, all)
    if (sum(ok) < max(20, ncol(Xfull_all) + 3)) next
    Xf <- Xfull_all[ok, , drop = FALSE]; Xp <- Xprot_all[ok, , drop = FALSE]
    fitf <- lm.fit(Xf, y[ok]); fitp <- lm.fit(Xp, y[ok])
    bf <- fitf$coefficients; bf[!is.finite(bf)] <- 0
    sse_full <- sum(fitf$residuals^2); sse_prot <- sum(fitp$residuals^2)
    pr2 <- if (sse_prot > 0) max(0, (sse_prot - sse_full) / sse_prot) else 0
    partial_r2[j] <- pr2
    use <- is.null(conditional_threshold) || pr2 >= conditional_threshold
    if (use && length(tech_cols)) {
      corrected[j, ok] <- y[ok] - as.numeric(Xf[, tech_cols, drop = FALSE] %*% bf[tech_cols])
      applied[j] <- TRUE
    }
  }
  list(matrix = corrected, partial_r2 = partial_r2, applied = applied)
}

partial_r2_term <- function(mat, meta, term) {
  ids <- intersect(colnames(mat), meta$sample_id)
  md <- meta[match(ids, meta$sample_id), , drop = FALSE]
  df <- data.frame(group2 = factor(md$group2), Run = factor(md$Run), Plate = factor(md$Plate))
  full <- model.matrix(~ group2 + Run + Plate, df)
  reduced <- model.matrix(switch(term, group2 = ~ Run + Plate, Run = ~ group2 + Plate, Plate = ~ group2 + Run), df)
  ans <- rep(NA_real_, nrow(mat))
  for (j in seq_len(nrow(mat))) {
    y <- mat[j, ids]; ok <- is.finite(y)
    if (sum(ok) < max(20, ncol(full) + 3)) next
    f1 <- lm.fit(full[ok, , drop = FALSE], y[ok]); f0 <- lm.fit(reduced[ok, , drop = FALSE], y[ok])
    s1 <- sum(f1$residuals^2); s0 <- sum(f0$residuals^2)
    ans[j] <- if (s0 > 0) max(0, (s0 - s1) / s0) else 0
  }
  ans
}

pca_summary <- function(mat, meta, max_features = 1000L) {
  ids <- intersect(colnames(mat), meta$sample_id); m <- mat[, ids, drop = FALSE]
  det <- rowMeans(is.finite(m)); v <- apply(m, 1, var, na.rm = TRUE)
  keep <- which(det >= 0.7 & is.finite(v) & v > 0)
  if (length(keep) > max_features) keep <- keep[order(v[keep], decreasing = TRUE)][seq_len(max_features)]
  x <- t(m[keep, , drop = FALSE])
  for (j in seq_len(ncol(x))) x[!is.finite(x[, j]), j] <- median(x[, j], na.rm = TRUE)
  pc <- prcomp(x, center = TRUE, scale. = TRUE)
  scores <- data.frame(sample_id = rownames(pc$x), PC1 = pc$x[, 1], PC2 = pc$x[, 2])
  md <- meta[match(scores$sample_id, meta$sample_id), , drop = FALSE]
  list(scores = cbind(scores, md[, c("group2", "Run", "Plate")]), variance = pc$sdev^2 / sum(pc$sdev^2), eta = c(PC1_group = eta2_factor(scores$PC1, md$group2), PC1_run = eta2_factor(scores$PC1, md$Run), PC1_plate = eta2_factor(scores$PC1, md$Plate), PC2_group = eta2_factor(scores$PC2, md$group2), PC2_run = eta2_factor(scores$PC2, md$Run), PC2_plate = eta2_factor(scores$PC2, md$Plate)))
}

pairwise_cor_median <- function(mat) {
  if (ncol(mat) < 2) return(NA_real_)
  cm <- suppressWarnings(cor(mat, use = "pairwise.complete.obs", method = "spearman"))
  median(cm[upper.tri(cm)], na.rm = TRUE)
}

control_metrics <- function(mat) c(median_feature_mad = median(apply(mat, 1, mad_or_na), na.rm = TRUE), sample_median_sd = sd(apply(mat, 2, median_or_na), na.rm = TRUE), median_pairwise_spearman = pairwise_cor_median(mat))
