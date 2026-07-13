plot_expr_boxplot <- function(mat,
                              metadata = NULL,
                              sample_id_col = "sample_id",
                              color_by = NULL,
                              log2_transform = FALSE,
                              add1 = 1,
                              drop_all_na_rows = FALSE,
                              show_outliers = FALSE,
                              show_points = FALSE,
                              show_median_label = FALSE,
                              median_label_digits = 2,
                              median_label_size = 3,
                              median_label_vjust = -0.5,
                              box_linewidth = 0.1,        # ★箱体边框+胡须线 线宽
                              median_linewidth = 0.5,      # ★中位数横线想保持的“绝对线宽”
                              point_alpha = 0.15,
                              point_size = 0.3,
                              order_by = c("data", "color_by", "none"),
                              facet_by = NULL,
                              ncol = NULL,
                              cols = NULL,
                              bold_title = FALSE,        # ★标题是否加粗
                              bold_axis_title = FALSE,   # ★坐标轴标题是否加粗
                              bold_legend_title = FALSE,  # ★图例标题是否加粗（可选）
                              fill_alpha = 0.7,          # ★颜色本身透明度（0~1）
                              box_alpha = 0.75,          # ★boxplot 几何层 alpha
                              panel_border_linewidth = 0.4,   # ★坐标轴边框线粗细
                              x_tick_linewidth = 0.2,    # ★横坐标刻度线粗细
                              y_tick_linewidth = 0.2,    # ★纵坐标刻度线粗细
                              title = NULL,
                              ylab = "Intensity",
                              font_family = "sans",
                              fs_title = 16,
                              fs_axis_title = 14,
                              fs_axis_text = 12,
                              fs_legend_title = 14,
                              fs_legend_text = 14,
                              show_x_text = TRUE,
                              show_x_ticks = TRUE,
                              ymin = NULL,
                              ymax = NULL,
                              xlab = "Sample",
                              transparent_bg = FALSE) {
  
  order_by <- match.arg(order_by)
  
  if (is.null(colnames(mat))) stop("mat 没有列名 colnames(mat)。")
  mat <- as.matrix(mat)
  
  if (isTRUE(drop_all_na_rows)) {
    keep_row <- rowSums(!is.na(mat)) > 0
    mat <- mat[keep_row, , drop = FALSE]
  }
  
  if (isTRUE(log2_transform)) {
    mat <- log2(mat + add1)
    if (is.null(ylab) || ylab == "Intensity") ylab <- paste0("log2(intensity + ", add1, ")")
  }
  
  df_long <- data.frame(mat, check.names = FALSE) |>
    tibble::rownames_to_column("Feature") |>
    tidyr::pivot_longer(cols = -Feature, names_to = "sample", values_to = "value") |>
    tidyr::drop_na(value)
  
  if (!is.null(metadata)) {
    if (!is.data.frame(metadata)) stop("metadata 必须是 data.frame。")
    if (!sample_id_col %in% colnames(metadata)) {
      stop(sprintf("metadata 中找不到 sample_id_col = '%s' 这一列。", sample_id_col))
    }
    
    need_cols <- unique(c(sample_id_col, color_by, facet_by))
    need_cols <- need_cols[!is.null(need_cols) & !is.na(need_cols)]
    meta_use <- metadata[, need_cols, drop = FALSE]
    meta_use[[sample_id_col]] <- as.character(meta_use[[sample_id_col]])
    
    df_long <- dplyr::left_join(
      df_long,
      meta_use,
      by = setNames(sample_id_col, "sample")
    )
    
    if (!is.null(color_by)) {
      miss2 <- unique(df_long$sample[is.na(df_long[[color_by]])])
      if (length(miss2) > 0) {
        warning("这些样本在 metadata 中未匹配到分组信息（将显示为 NA）：\n",
                paste(miss2, collapse = ", "))
      }
    }
  }
  
  # ---------- 样本顺序 ----------
  grp_levels <- NULL   # ★新增：用于控制图例顺序
  sample_levels <- colnames(mat)
  
  if (!is.null(metadata) && !is.null(color_by) && order_by == "color_by") {
    
    meta_map <- metadata[, c(sample_id_col, color_by), drop = FALSE]
    meta_map[[sample_id_col]] <- as.character(meta_map[[sample_id_col]])
    meta_map[[color_by]]      <- as.character(meta_map[[color_by]])
    
    ord_df <- data.frame(
      sample = colnames(mat),
      orig   = seq_along(colnames(mat)),
      stringsAsFactors = FALSE
    )
    
    idx <- match(ord_df$sample, meta_map[[sample_id_col]])
    ord_df$grp <- meta_map[[color_by]][idx]
    
    # 没匹配到的放最后
    if (any(is.na(ord_df$grp))) {
      warning("这些样本在 metadata 中没有匹配到分组信息，将排在最后：\n",
              paste(ord_df$sample[is.na(ord_df$grp)], collapse = ", "))
      ord_df$grp[is.na(ord_df$grp)] <- "NA"
    }
    
    # 如果当前顺序已经“按组连续”，就不再重排
    grp_seq <- ord_df$grp
    r <- rle(grp_seq)
    already_grouped <- !any(duplicated(r$values))
    
    if (!already_grouped) {
      # 需要排序：按组排序，组内保持原顺序
      ord_df <- ord_df[order(ord_df$grp, ord_df$orig), , drop = FALSE]
    }
    
    sample_levels <- ord_df$sample
  }
  
  df_long$sample <- factor(df_long$sample, levels = sample_levels)
  
  # ---------- 构造 ggplot ----------
  legend_labels <- NULL  # ★新增：图例显示 n 的标签
  
  if (!is.null(metadata) && !is.null(color_by)) {
    
    # ★把 NA 分组显式变成 "NA"，这样也能计数并显示在图例里
    # 先拿到 metadata 里 Plate 的 levels（如果它是 factor）
    lv <- NULL
    if (is.factor(metadata[[color_by]])) {
      lv <- levels(metadata[[color_by]])
    }
    
    df_long[[color_by]] <- as.character(df_long[[color_by]])
    if (anyNA(df_long[[color_by]])) df_long[[color_by]][is.na(df_long[[color_by]])] <- "NA"
    
    # 用 metadata 的 levels 固定顺序
    if (!is.null(lv)) {
      if (!("NA" %in% lv) && any(df_long[[color_by]] == "NA")) lv <- c(lv, "NA")
      df_long[[color_by]] <- factor(df_long[[color_by]], levels = lv)
    } else {
      df_long[[color_by]] <- factor(df_long[[color_by]])
    }
    
    
    # ★计算每组样本量（按 sample 去重）
    n_tbl <- df_long |>
      dplyr::distinct(sample, .data[[color_by]]) |>
      dplyr::count(.data[[color_by]], name = "n") |>
      dplyr::mutate(
        .grp = as.character(.data[[color_by]]),
        .lbl = paste0(.grp, " (n=", n, ")")
      )
    
    legend_labels <- stats::setNames(n_tbl$.lbl, n_tbl$.grp)
    
    p <- ggplot2::ggplot(df_long, ggplot2::aes(x = sample, y = value, fill = .data[[color_by]]))
    fill_levels <- levels(df_long[[color_by]])
    
    # 补齐可能存在但没出现的 level（n=0）
    miss_lv <- setdiff(fill_levels, names(legend_labels))
    if (length(miss_lv) > 0) {
      legend_labels[miss_lv] <- paste0(miss_lv, " (n=0)")
    }
    legend_labels <- legend_labels[fill_levels]
    
  } else {
    p <- ggplot2::ggplot(df_long, ggplot2::aes(x = sample, y = value, fill = sample))
    fill_levels <- levels(df_long$sample)
  }
  
  # ---------- boxplot ----------
  if (!is.numeric(box_linewidth) || box_linewidth <= 0) stop("box_linewidth 必须是 >0 的数值。")
  if (!is.numeric(median_linewidth) || median_linewidth <= 0) stop("median_linewidth 必须是 >0 的数值。")
  
  # 让中位数横线保持固定粗细：median_linewidth = box_linewidth * fatten
  fatten_use <- median_linewidth / box_linewidth
  
  p <- p + ggplot2::geom_boxplot(
    width = 0.6,
    outlier.shape = if (isTRUE(show_outliers)) 16 else NA,
    alpha = box_alpha,
    linewidth = box_linewidth,
    fatten = fatten_use
  )
  
  if (isTRUE(show_median_label)) {
    p <- p + ggplot2::stat_summary(
      fun = median,
      geom = "text",
      ggplot2::aes(
        label = sprintf(
          paste0("%.", median_label_digits, "f"),
          after_stat(y)
        )
      ),
      vjust = median_label_vjust,
      size = median_label_size,
      color = "black",
      show.legend = FALSE
    )
  }
  
  if (isTRUE(show_points)) {
    p <- p + ggplot2::geom_jitter(width = 0.18, height = 0,
                                  alpha = point_alpha, size = point_size, color = "black")
  }
  
  if (!is.null(metadata) && !is.null(facet_by)) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste0("~", facet_by)),
                                 scales = "free_x", ncol = ncol)
  }
  
  # ---------- 颜色 + “颜色本身透明度(fill_alpha)” ----------
  apply_alpha <- function(col_vec, a) {
    a <- as.numeric(a)
    if (!is.finite(a) || a < 0 || a > 1) stop("fill_alpha 必须是 0~1 之间的数值。")
    if (requireNamespace("scales", quietly = TRUE)) {
      return(scales::alpha(col_vec, a))
    } else {
      return(grDevices::adjustcolor(col_vec, alpha.f = a))
    }
  }
  
  # 组装 colors：优先用用户 cols；否则在需要透明度时自动生成一套
  need_manual_scale <- !is.null(cols) || (!isTRUE(all.equal(fill_alpha, 1)))
  
  if (need_manual_scale) {
    
    if (is.null(cols)) {
      # 自动生成颜色（为了能应用 fill_alpha）
      if (!requireNamespace("scales", quietly = TRUE)) {
        cols_use <- grDevices::hcl.colors(length(fill_levels), palette = "Dynamic")
      } else {
        cols_use <- scales::hue_pal()(length(fill_levels))
      }
      cols_use <- stats::setNames(cols_use, fill_levels)
      
    } else {
      cols <- as.character(cols)
      has_names <- !is.null(names(cols)) && any(nzchar(names(cols)))
      
      if (has_names) {
        miss <- setdiff(fill_levels, names(cols))
        if (length(miss) > 0) {
          stop("cols 缺少这些水平的颜色：", paste(miss, collapse = ", "),
               "\n当前 fill_levels: ", paste(fill_levels, collapse = ", "))
        }
        cols_use <- cols[fill_levels]
      } else {
        if (length(cols) < length(fill_levels)) {
          stop("cols 颜色数量不足：给了 ", length(cols), " 个，但需要 ", length(fill_levels), " 个。")
        }
        cols_use <- cols[seq_along(fill_levels)]
        names(cols_use) <- fill_levels
      }
    }
    
    # 颜色本身加透明度
    cols_use <- apply_alpha(cols_use, fill_alpha)
    
    p <- p + ggplot2::scale_fill_manual(
      values = cols_use,
      drop   = FALSE,
      labels = if (!is.null(legend_labels)) legend_labels else ggplot2::waiver()
    )
    
  } else {
    
    # 不需要 manual scale：保留你原来的逻辑，并在有分组时加 labels（显示 n）
    if (!is.null(metadata) && !is.null(color_by) && !is.null(legend_labels)) {
      p <- p + ggplot2::scale_fill_discrete(
        drop   = FALSE,
        labels = legend_labels
      )
    } else {
      if (is.null(metadata) || is.null(color_by)) {
        p <- p + ggplot2::guides(fill = "none")
      }
    }
  }
  
  # ---------- 背景 ----------
  bg_fill <- if (isTRUE(transparent_bg)) "transparent" else "white"
  
  p <- p +
    ggplot2::labs(title = title, x = xlab, y = ylab, fill = color_by) +
    ggplot2::theme_minimal(base_family = font_family) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      
      plot.title  = ggplot2::element_text(size = fs_title, hjust = 0.5,
                                          face = if (bold_title) "bold" else "plain"),
      
      axis.title.x = ggplot2::element_text(size = fs_axis_title,
                                           face = if (bold_axis_title) "bold" else "plain"),
      axis.title.y = ggplot2::element_text(size = fs_axis_title,
                                           face = if (bold_axis_title) "bold" else "plain"),
      
      axis.text.x  = if (show_x_text)
        ggplot2::element_text(size = fs_axis_text, angle = 45, hjust = 1, color = "black")
      else
        ggplot2::element_blank(),
      
      axis.ticks.x = if (show_x_ticks) ggplot2::element_line(linewidth = x_tick_linewidth) else ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_line(linewidth = y_tick_linewidth),
      
      axis.text.y = ggplot2::element_text(size = fs_axis_text, color = "black"),
      
      legend.title = ggplot2::element_text(size = fs_legend_title,
                                           face = if (bold_legend_title) "bold" else "plain"),
      legend.text  = ggplot2::element_text(size = fs_legend_text),
      
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = panel_border_linewidth),
      
      panel.background = ggplot2::element_rect(fill = bg_fill, colour = NA),
      plot.background  = ggplot2::element_rect(fill = bg_fill, colour = NA),
      legend.background = ggplot2::element_rect(fill = bg_fill, colour = NA),
      legend.box.background = ggplot2::element_rect(fill = bg_fill, colour = NA)
    )
  
  if (!is.null(ymin) || !is.null(ymax)) {
    p <- p + ggplot2::coord_cartesian(
      ylim = c(ifelse(is.null(ymin), NA, ymin),
               ifelse(is.null(ymax), NA, ymax))
    )
  }
  
  return(p)
}
