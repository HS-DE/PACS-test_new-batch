align_samples <- function(exp_data,
                          group_data,
                          sample_col   = "sample_id",
                          group_col    = "group",          # 可以是一个字符，也可以是 c("Study","group")
                          keep_groups  = NULL,             # 单列：向量；多列：named list
                          group_levels = NULL,             # 单列：向量；多列：named list
                          sort_within  = c("metadata", "sample_name"),
                          verbose      = TRUE) {
  
  sort_within <- match.arg(sort_within)
  group_col <- as.character(group_col)
  
  # -----------------------------
  # 0) 基础检查
  # -----------------------------
  need_cols <- unique(c(sample_col, group_col))
  if (!all(need_cols %in% names(group_data))) {
    stop("group_data 中必须包含列: ", paste(need_cols, collapse = ", "))
  }
  if (is.null(colnames(exp_data))) {
    stop("exp_data 必须有列名（列名=样本ID）。")
  }
  
  # -----------------------------
  # 1) 取必要列，转字符；样本去重（保留第一次）
  # -----------------------------
  meta <- group_data[, need_cols, drop = FALSE]
  meta[[sample_col]] <- as.character(meta[[sample_col]])
  for (gc in group_col) meta[[gc]] <- as.character(meta[[gc]])
  
  if (anyDuplicated(meta[[sample_col]]) > 0) {
    if (verbose) message("发现 group_data 中样本ID重复：保留首次出现，丢弃其余。")
    meta <- meta[!duplicated(meta[[sample_col]]), , drop = FALSE]
  }
  
  # -----------------------------
  # 2) 两边样本交集，只保留交集
  # -----------------------------
  samples_intersect <- intersect(meta[[sample_col]], colnames(exp_data))
  dropped_from_meta <- setdiff(meta[[sample_col]], colnames(exp_data))
  dropped_from_exp  <- setdiff(colnames(exp_data), meta[[sample_col]])
  
  if (length(samples_intersect) == 0) stop("样本交集为空：请检查样本ID是否一致。")
  meta <- meta[meta[[sample_col]] %in% samples_intersect, , drop = FALSE]
  
  # -----------------------------
  # 3) 过滤 keep_groups
  #    单列：keep_groups 是向量
  #    多列：keep_groups 是 named list，如 list(Study=c("S1","S2"), group=c("HC","CRC"))
  # -----------------------------
  if (!is.null(keep_groups)) {
    if (length(group_col) == 1 && !is.list(keep_groups)) {
      meta <- meta[meta[[group_col]] %in% keep_groups, , drop = FALSE]
    } else {
      if (!is.list(keep_groups)) stop("当 group_col 有多个时，keep_groups 必须是 named list。")
      for (gc in group_col) {
        if (!is.null(keep_groups[[gc]])) {
          meta <- meta[meta[[gc]] %in% keep_groups[[gc]], , drop = FALSE]
        }
      }
    }
    if (nrow(meta) == 0) stop("过滤 keep_groups 后没有样本。")
  }
  
  # -----------------------------
  # 4) 为每个 group_col 设定 levels（用于排序）
  #    单列：group_levels 向量
  #    多列：group_levels named list
  # -----------------------------
  if (length(group_col) == 1 && !is.null(group_levels) && !is.list(group_levels)) {
    # 单列、用户给的是向量
    gc <- group_col
    meta[[gc]] <- factor(meta[[gc]], levels = group_levels)
  } else {
    # 多列（或用户传了 list）
    if (is.null(group_levels)) {
      # 没给 levels：默认按 meta 中出现顺序
      group_levels <- lapply(group_col, function(gc) unique(meta[[gc]]))
      names(group_levels) <- group_col
    } else {
      if (!is.list(group_levels)) stop("当 group_col 有多个时，group_levels 必须是 named list。")
      # 没指定的列，默认按出现顺序
      for (gc in group_col) {
        if (is.null(group_levels[[gc]])) group_levels[[gc]] <- unique(meta[[gc]])
      }
    }
    for (gc in group_col) {
      meta[[gc]] <- factor(meta[[gc]], levels = group_levels[[gc]])
    }
  }
  
  # 去掉不在 levels 中的（会变 NA）
  for (gc in group_col) meta <- meta[!is.na(meta[[gc]]), , drop = FALSE]
  if (nrow(meta) == 0) stop("所有样本的组别都不在 group_levels 中。")
  
  # -----------------------------
  # 5) 排序：支持多级排序
  #    先按 group_col[1]（如 Study），再按 group_col[2]（如 group）...
  #    组内：metadata 顺序（稳定）或 sample_name
  # -----------------------------
  meta$.orig <- seq_len(nrow(meta))  # 用来保持组内稳定
  
  order_keys <- lapply(group_col, function(gc) meta[[gc]])
  if (sort_within == "sample_name") {
    order_keys <- c(order_keys, list(meta[[sample_col]]))
  } else {
    order_keys <- c(order_keys, list(meta$.orig))
  }
  ord <- do.call(order, order_keys)
  meta <- meta[ord, , drop = FALSE]
  meta$.orig <- NULL
  
  # -----------------------------
  # 6) 重排表达矩阵列
  # -----------------------------
  exp_aligned <- exp_data[, meta[[sample_col]], drop = FALSE]
  stopifnot(identical(colnames(exp_aligned), meta[[sample_col]]))
  
  if (verbose) {
    message(sprintf("样本对齐完成：%d 个样本。", ncol(exp_aligned)))
    if (length(dropped_from_meta) > 0)
      message(sprintf("  meta里有但表达矩阵缺失：%d 个", length(dropped_from_meta)))
    if (length(dropped_from_exp) > 0)
      message(sprintf("  表达矩阵里有但meta缺失：%d 个", length(dropped_from_exp)))
  }
  
  return(list(
    exp               = exp_aligned,
    meta              = meta,
    dropped_from_meta = dropped_from_meta,
    dropped_from_exp  = dropped_from_exp
  ))
}

if(F){
  out <- align_samples(
    exp_data     = matrix2,
    group_data   = metadata3,
    sample_col   = "sample_id",
    group_col    = c("Study", "group"),
    keep_groups  = list(Study = c("S1","S2")),
    group_levels = list(
      Study = c("S1","S2"),
      group = c("HC","CRC")   # 你想要的组顺序
    ),
    sort_within  = "metadata",
    verbose      = TRUE
  )
  
  
}