# 删除 fviz_pca_ind 默认加的 x=0/y=0 参考线（如果存在）
remove_origin_lines <- function(p) {
  keep <- vapply(p$layers, function(l) {
    !(inherits(l$geom, "GeomHline") || inherits(l$geom, "GeomVline"))
  }, logical(1))
  p$layers <- p$layers[keep]
  p
}

# 用“和你 plot_pca_samples 同样的 PCA 逻辑”计算多个矩阵的统一坐标范围
get_common_pca_limits <- function(mat_list, axes = c(1,2),
                                  symmetric = TRUE, pad_frac = 0.05) {
  coords_all <- do.call(rbind, lapply(mat_list, function(mat_expr) {
    
    # 1) 转成 样本 x 特征
    X <- t(as.matrix(mat_expr))
    X[!is.finite(X)] <- NA
    
    # 2) 过滤特征：至少 2 个非 NA + 方差>0
    keep1 <- colSums(!is.na(X)) >= 2
    X <- X[, keep1, drop = FALSE]
    v <- apply(X, 2, var, na.rm = TRUE)
    keep2 <- is.finite(v) & v > 0
    X <- X[, keep2, drop = FALSE]
    
    # 3) PCA（与主函数一致）
    pca_res <- FactoMineR::PCA(as.data.frame(X), scale.unit = TRUE,
                               ncp = max(5, max(axes)), graph = FALSE)
    
    pca_res$ind$coord[, axes, drop = FALSE]
  }))
  
  if (symmetric) {
    m <- max(abs(coords_all), na.rm = TRUE)
    m <- m * (1 + pad_frac)
    list(xlim = c(-m, m), ylim = c(-m, m))
  } else {
    pad <- function(r) { d <- diff(r); r + c(-pad_frac*d, pad_frac*d) }
    list(xlim = pad(range(coords_all[,1], na.rm=TRUE)),
         ylim = pad(range(coords_all[,2], na.rm=TRUE)))
  }
}


plot_pca_samples <- function(mat_expr,
                             group,
                             group_levels = NULL,
                             metadata = NULL,                 # ★新增：可选 metadata
                             sample_col = "sample_id",         # ★新增：metadata 里样本ID列名
                             title = NULL,
                             pointsize = 3,
                             axes = c(1, 2),
                             addEllipses = TRUE,
                             ellipse.level = 0.9,
                             palette = "aaas",          # "aaas" / "hue" / 颜色向量(命名或不命名)
                             palette_alpha = 1,         # ★颜色本身透明度（0~1）
                             mean.point = FALSE,
                             fs_title = 16,
                             fs_legend = 14,
                             fs_axis = 12,
                             bold_title  = FALSE,
                             bold_axis   = FALSE,
                             bold_legend = FALSE,
                             show_legend = TRUE,
                             legend_position = c("right","left","top","bottom","none"),
                             
                             show_labels = FALSE,
                             label_size = 3,
                             transparent_bg = FALSE,
                             max_overlaps = 50,
                             axis_title_style = c("PC", "blank"),  # PC: 显示PC1/PC2；blank: 不显示
                             xlim = NULL,
                             ylim = NULL,
                             axis_on_edge = TRUE,                 # 轴线画在边缘
                             remove_origin = TRUE,                # 去掉x=0/y=0参考线
                             equal_aspect = FALSE                 # TRUE时 x/y 同比例显示（可选）
) {
  
  # =========================
  # 0) 依赖检查
  # =========================
  need_pkgs <- c("FactoMineR", "factoextra", "ggplot2")
  miss <- need_pkgs[!vapply(need_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss) > 0) {
    stop("缺少依赖包：", paste(miss, collapse = ", "),
         "\n请先安装：install.packages(c(", paste0("'", miss, "'", collapse = ", "), "))")
  }
  if (isTRUE(show_labels) && !requireNamespace("ggrepel", quietly = TRUE)) {
    stop("show_labels=TRUE 需要 ggrepel：install.packages('ggrepel')")
  }
  
  # palette_alpha 检查
  palette_alpha <- as.numeric(palette_alpha)
  if (!is.finite(palette_alpha) || palette_alpha < 0 || palette_alpha > 1) {
    stop("palette_alpha 必须是 0~1 之间的数值。")
  }
  
  if (missing(group_levels)) {
    group_levels <- get0("group_levels", envir = parent.frame(), inherits = TRUE, ifnotfound = NULL)
  }
  if (!is.null(group_levels)) {
    group_levels <- as.character(group_levels)
    group_levels <- group_levels[!is.na(group_levels) & nzchar(group_levels)]
    group_levels <- group_levels[!duplicated(group_levels)]
    if (length(group_levels) == 0) group_levels <- NULL
  }
  
  # =========================
  # 0.5) 小工具：去掉原点参考线（去掉所有 hline/vline 层）
  # =========================
  remove_origin_lines <- function(p) {
    if (length(p$layers) == 0) return(p)
    keep <- vapply(p$layers, function(l) {
      geom_class <- class(l$geom)[1]
      !(geom_class %in% c("GeomHline", "GeomVline"))
    }, logical(1))
    p$layers <- p$layers[keep]
    p
  }
  
  # =========================
  # 1) 转成：样本 x 特征（输入通常是 特征 x 样本，所以这里 t()）
  # =========================
  orig_sample_ids <- colnames(mat_expr)
  if (is.null(orig_sample_ids) || any(!nzchar(orig_sample_ids))) {
    stop("mat_expr 必须有列名（样本ID）。请先给 mat_expr 设置 colnames（样本名）。")
  }
  
  mat_expr <- as.data.frame(t(mat_expr), check.names = FALSE, stringsAsFactors = FALSE)
  
  # 把样本名放到 rownames（用于后续可视化/标签）
  if (!is.null(orig_sample_ids) && length(orig_sample_ids) == nrow(mat_expr)) {
    rownames(mat_expr) <- as.character(orig_sample_ids)
  }
  
  n_samples <- nrow(mat_expr)
  
  # =========================
  # 1.5) ★新增：对齐 group 到 mat_expr 的样本顺序
  # 支持两种方式：
  #   A) group = "Study"（metadata 列名）
  #   B) group = 向量（可无names），配合 metadata + sample_col 对齐
  # =========================
  align_group_by_metadata <- function(group, metadata, sample_col, sample_ids) {
    # group 是列名
    if (is.character(group) && length(group) == 1 && group %in% colnames(metadata)) {
      gv <- metadata[[group]]
      names(gv) <- sample_ids
      return(gv)
    }
    
    # group 是向量
    gv <- group
    if (is.null(names(gv)) || all(!nzchar(names(gv)))) {
      # 如果长度等于 nrow(metadata)，则假设 gv 与 metadata 行对应
      if (length(gv) == nrow(metadata)) {
        names(gv) <- sample_ids
      }
    }
    return(gv)
  }
  
  if (!is.null(metadata)) {
    metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)
    
    if (!(sample_col %in% colnames(metadata))) {
      stop("metadata 中找不到 sample_col = '", sample_col, "' 这一列。\n",
           "metadata 的列名包括：", paste(colnames(metadata), collapse = ", "))
    }
    
    sample_ids_meta <- as.character(metadata[[sample_col]])
    if (any(!nzchar(sample_ids_meta))) {
      stop("metadata$", sample_col, " 存在空字符串样本ID，请先清理。")
    }
    if (anyDuplicated(sample_ids_meta)) {
      dup_ids <- unique(sample_ids_meta[duplicated(sample_ids_meta)])
      stop("metadata$", sample_col, " 存在重复样本ID：", paste(dup_ids, collapse = ", "),
           "\n请先去重（一个样本ID只能对应一条分组信息）。")
    }
    
    # 先把 group 规范成“带 names 的向量”，names=样本ID
    group_named <- align_group_by_metadata(group, metadata, sample_col, sample_ids_meta)
    
    # 如果现在仍然没有 names，就无法按样本名对齐
    if (is.null(names(group_named)) || all(!nzchar(names(group_named)))) {
      stop("你传入了 metadata，但 group 既不是 metadata 的列名，也不是可对齐的向量。\n",
           "推荐用法：group='分组列名'，例如 group='Study'；或 group=metadata$Study。")
    }
    
    # 按 mat_expr 的样本顺序对齐
    group_aligned <- group_named[orig_sample_ids]
    miss_ids <- orig_sample_ids[is.na(group_aligned)]
    if (length(miss_ids) > 0) {
      stop("以下样本在 metadata$", sample_col, " 中找不到对应分组信息：\n",
           paste(miss_ids, collapse = ", "), "\n",
           "请检查：mat_expr 的列名是否与 metadata$", sample_col, " 完全一致。")
    }
    
    group <- group_aligned
  } else {
    # 没给 metadata：如果 group 是“带names向量”，也尝试按样本名对齐（可选但很有用）
    if (!is.null(names(group)) && all(orig_sample_ids %in% names(group))) {
      group <- group[orig_sample_ids]
    }
  }
  
  # 对齐后再做长度检查（现在就不会“长度对了但顺序错了”）
  if (length(group) != n_samples) {
    stop("group 长度 (", length(group), ") 必须等于样本数 (", n_samples, ")。\n",
         "提示：mat_expr 应为 特征x样本；函数内部会 t() 成 样本x特征。\n",
         "如果你使用 metadata 对齐，请确认 mat_expr 的列名是样本ID。")
  }
  
  # =========================
  # 2) 清理非有限值 -> NA
  # =========================
  X <- as.matrix(mat_expr)
  X[!is.finite(X)] <- NA
  
  # =========================
  # 3) 过滤特征：
  #    - 至少 2 个非 NA
  #    - 方差 > 0
  # =========================
  keep1 <- colSums(!is.na(X)) >= 2
  X <- X[, keep1, drop = FALSE]
  if (ncol(X) < 2) stop("过滤后特征数 < 2，无法做 PCA（可能 NA 太多）。")
  
  v <- apply(X, 2, var, na.rm = TRUE)
  keep2 <- is.finite(v) & v > 0
  X <- X[, keep2, drop = FALSE]
  if (ncol(X) < 2) stop("过滤后方差>0的特征数 < 2，无法做 PCA。")
  
  # =========================
  # 4) PCA（FactoMineR 会对 NA 做均值插补）
  # =========================
  max_axis <- max(axes)
  pca_res <- FactoMineR::PCA(as.data.frame(X), scale.unit = TRUE, ncp = max(5, max_axis), graph = FALSE)
  
  # =========================
  # 5) 分组（factor）+ 计算每组 n，并拼到图例标签
  #    ★关键：保留用户传入 group 的 levels
  # =========================
  orig_levels <- if (is.factor(group)) levels(group) else NULL
  
  group_chr <- as.character(group)
  group_chr[is.na(group_chr)] <- "NA"   # NA 也作为一组
  
  if (!is.null(group_levels)) {
    levs <- group_levels
    missing_levs <- setdiff(unique(group_chr[group_chr != "NA"]), levs)
    if (length(missing_levs) > 0) levs <- c(levs, missing_levs)
    if (any(group_chr == "NA") && !("NA" %in% levs)) levs <- c(levs, "NA")
    group_fac <- factor(group_chr, levels = levs)
  } else if (!is.null(orig_levels)) {
    levs <- orig_levels
    if (any(group_chr == "NA") && !("NA" %in% levs)) levs <- c(levs, "NA")
    group_fac <- factor(group_chr, levels = levs)
  } else {
    levs <- unique(group_chr)
    levs_no_na <- levs[levs != "NA"]
    num_try <- suppressWarnings(as.numeric(levs_no_na))
    if (length(levs_no_na) > 0 && all(!is.na(num_try))) {
      levs_no_na <- as.character(sort(num_try))
    } else {
      levs_no_na <- sort(levs_no_na)
    }
    if ("NA" %in% levs) levs <- c(levs_no_na, "NA") else levs <- levs_no_na
    group_fac <- factor(group_chr, levels = levs)
  }
  
  group_levels_used <- levels(group_fac)
  n_groups <- length(group_levels_used)
  
  n_vec <- as.integer(table(group_fac)[group_levels_used])
  legend_labels <- stats::setNames(
    paste0(group_levels_used, " (n=", n_vec, ")"),
    group_levels_used
  )
  
  # =========================
  # 6) 先画基础 PCA 图
  # =========================
  p <- factoextra::fviz_pca_ind(
    pca_res,
    pointsize = pointsize,
    axes = axes,
    label = "none",
    addEllipses = addEllipses,
    ellipse.level = ellipse.level,
    habillage = group_fac,
    mean.point = mean.point,
    title = title
  )
  
  # =========================
  # 7) 只改“颜色本身透明度”：给颜色向量加 alpha，再用 manual scale
  # =========================
  add_color_alpha <- function(col_vec, a) {
    if (isTRUE(all.equal(a, 1))) return(col_vec)
    col_vec <- as.character(col_vec)
    if (requireNamespace("scales", quietly = TRUE)) {
      return(scales::alpha(col_vec, a))
    } else {
      return(grDevices::adjustcolor(col_vec, alpha.f = a))
    }
  }
  
  get_pal_use <- function(palette, levs, n_groups) {
    if (is.character(palette) && length(palette) > 1) {
      pal_vec <- as.character(palette)
      has_names <- !is.null(names(pal_vec)) && any(nzchar(names(pal_vec)))
      if (!has_names) {
        if (length(pal_vec) < length(levs)) {
          stop("你提供的颜色数量 (", length(pal_vec),
               ") 少于组数 (", length(levs), ")。\n",
               "请补足颜色，或使用命名向量：c(group1='#..', group2='#..', ...)")
        }
        pal_use <- pal_vec[seq_along(levs)]
        names(pal_use) <- levs
      } else {
        missing <- setdiff(levs, names(pal_vec))
        if (length(missing) > 0) {
          stop("palette 命名向量缺少这些组的颜色：", paste(missing, collapse = ", "), "\n",
               "你的组 levels: ", paste(levs, collapse = ", "))
        }
        pal_use <- pal_vec[levs]
      }
      return(pal_use)
    }
    
    if (!(is.character(palette) && length(palette) == 1)) {
      stop("palette 只能是：\n",
           "1) 单个字符串：'aaas' 或 'hue'\n",
           "2) 颜色向量：c('#E64B35', '#4DBBD5', ...)\n",
           "3) 命名颜色向量：c(A='#E64B35', B='#4DBBD5', ...)")
    }
    
    pal_str <- tolower(palette)
    
    if (pal_str == "aaas" && requireNamespace("ggsci", quietly = TRUE) && n_groups <= 10) {
      pal_use <- ggsci::pal_aaas("default")(n_groups)
      names(pal_use) <- levs
      return(pal_use)
    }
    
    if (requireNamespace("scales", quietly = TRUE)) {
      pal_use <- scales::hue_pal()(n_groups)
    } else {
      pal_use <- grDevices::hcl.colors(n_groups, palette = "Dynamic")
    }
    names(pal_use) <- levs
    return(pal_use)
  }
  
  pal_use <- get_pal_use(palette, group_levels_used, n_groups)
  pal_use <- add_color_alpha(pal_use, palette_alpha)
  
  legend_title <- NULL
  
  p <- p +
    ggplot2::scale_color_manual(
      values = pal_use,
      breaks = group_levels_used,
      labels = legend_labels[group_levels_used],
      drop = FALSE,
      name = legend_title
    ) +
    ggplot2::scale_fill_manual(
      values = pal_use,
      breaks = group_levels_used,
      labels = legend_labels[group_levels_used],
      drop = FALSE,
      name = legend_title
    )
  
  p <- p + ggplot2::scale_shape_manual(
    values = stats::setNames(rep(16, length(group_levels_used)), group_levels_used),
    breaks = group_levels_used,
    labels = legend_labels[group_levels_used],
    drop = FALSE,
    guide = "none",
    name = legend_title
  )
  
  # =========================
  # 8) 背景与主题
  # =========================
  bg_fill <- if (isTRUE(transparent_bg)) "transparent" else "white"
  
  p <- p +
    ggplot2::labs(title = title) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      
      panel.background = ggplot2::element_rect(fill = bg_fill, colour = NA),
      plot.background  = ggplot2::element_rect(fill = bg_fill, colour = NA),
      legend.background = ggplot2::element_rect(fill = bg_fill, colour = NA),
      legend.box.background = ggplot2::element_rect(fill = bg_fill, colour = NA),
      
      plot.title = ggplot2::element_text(
        size = fs_title,
        face = if (isTRUE(bold_title)) "bold" else "plain"
      ),
      legend.title = ggplot2::element_text(
        size = fs_legend,
        face = if (isTRUE(bold_legend)) "bold" else "plain"
      ),
      legend.text = ggplot2::element_text(
        size = fs_legend,
        face = if (isTRUE(bold_legend)) "bold" else "plain"
      ),
      axis.text = ggplot2::element_text(
        size = fs_axis,
        face = if (isTRUE(bold_axis)) "bold" else "plain"
      ),
      axis.title = ggplot2::element_text(
        size = fs_axis,
        face = if (isTRUE(bold_axis)) "bold" else "plain"
      )
    )
  
  # =========================
  # 8.5) 坐标轴标题：Dim -> PC（或直接隐藏）
  # =========================
  axis_title_style <- match.arg(axis_title_style)
  pc_pct <- round(pca_res$eig[axes, 2], 1)
  
  if (axis_title_style == "PC") {
    xlab <- sprintf("PC%d (%.1f%%)", axes[1], pc_pct[1])
    ylab2 <- sprintf("PC%d (%.1f%%)", axes[2], pc_pct[2])
  } else {
    xlab <- NULL; ylab2 <- NULL
  }
  p <- p + ggplot2::labs(x = xlab, y = ylab2)
  
  # =========================
  # 8.6) 固定坐标范围
  # =========================
  if (!is.null(xlim) || !is.null(ylim)) {
    if (isTRUE(equal_aspect)) {
      p <- p + ggplot2::coord_equal(xlim = xlim, ylim = ylim, expand = FALSE)
    } else {
      p <- p + ggplot2::coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE)
    }
  }
  
  # =========================
  # 8.7) 轴线在边缘 + 去掉原点参考线
  # =========================
  if (isTRUE(remove_origin)) {
    p <- remove_origin_lines(p)
  }
  if (isTRUE(axis_on_edge)) {
    p <- p + ggplot2::theme(
      axis.line  = ggplot2::element_line(),
      axis.ticks = ggplot2::element_line()
    )
  }
  
  # =========================
  # 8.8) 图例开关
  # =========================
  legend_position <- match.arg(legend_position)
  
  if (!isTRUE(show_legend) || legend_position == "none") {
    p <- p + ggplot2::theme(legend.position = "none")
  } else {
    p <- p + ggplot2::theme(legend.position = legend_position)
  }
  
  # =========================
  # 9) 标签
  # =========================
  if (isTRUE(show_labels)) {
    p <- p + ggrepel::geom_text_repel(
      ggplot2::aes(label = name),
      show.legend = FALSE,
      size = label_size,
      max.overlaps = max_overlaps
    )
  }
  
  return(p)
}
if(F){
  p <- plot_pca_samples(
    mat_expr,
    group = "Study",
    metadata = metadata3,
    sample_col = "sample_id",
    title = "PCA"
  )
  print(p)
}
