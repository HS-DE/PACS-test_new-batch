plot_cv_by_metadata <- function(
  mat,
  metadata,
  sample_id_col = "sample_id",
  group_col = "group",
  group_levels = NULL,
  cols = NULL,
  na_rm = TRUE,
  min_non_na = 3,
  min_detect_rate = 0.70,
  min_mean = 1e-8,
  cv_cap = NULL,
  add_violin = TRUE,
  add_box = TRUE,
  angle = NULL,
  hjust = 0.5,
  vjust = NULL,

  # ===== 新增：散点参数 =====
  add_points = FALSE,
  point_alpha = 0.8,
  point_size = 1.8,
  point_width = 0.15,
  point_border = TRUE,
  point_stroke = 0.1,

  title = NULL,
  subtitle = NULL,
  ylab = "Coefficient of Variation (%)",
  font_family = NULL,
  fs_title = 16,
  fs_subtitle = 12,
  fs_axis_title = 14,
  fs_axis_text = 12,
  fs_strip_text = 20,
  panel_border_size = 1,
  transparent_bg = FALSE,

  bold_title = FALSE,
  bold_axis_title = FALSE,
  bold_strip_text = FALSE,
  bold_median_label = FALSE,

  out_prefix = NULL,
  save_formats = NULL,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300,
  tiff_compression = "lzw"
) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  # =========================
  # 0) 输入检查
  # =========================
  if (is.null(colnames(mat))) {
    stop("mat 没有列名 colnames(mat)。")
  }
  mat <- as.matrix(mat)

  if (!is.data.frame(metadata)) {
    stop("metadata 必须是 data.frame。")
  }
  if (!sample_id_col %in% colnames(metadata)) {
    stop(sprintf("metadata 中找不到 sample_id_col='%s'。", sample_id_col))
  }
  if (!group_col %in% colnames(metadata)) {
    stop(sprintf("metadata 中找不到 group_col='%s'。", group_col))
  }

  # =========================
  # 1) 建立 sample -> group 映射，并过滤到 mat 中存在的样本
  # =========================
  map_df <- metadata[, c(sample_id_col, group_col), drop = FALSE]
  map_df[[sample_id_col]] <- as.character(map_df[[sample_id_col]])
  map_df[[group_col]] <- as.character(map_df[[group_col]])

  # 只保留 mat 里存在的列
  map_df <- map_df[map_df[[sample_id_col]] %in% colnames(mat), , drop = FALSE]

  miss_samples <- setdiff(colnames(mat), map_df[[sample_id_col]])
  if (length(miss_samples) > 0) {
    warning(
      "这些样本在 metadata 中找不到，将不会参与 CV 计算：\n",
      paste(miss_samples, collapse = ", ")
    )
  }

  # 每组样本列名列表
  group_cols <- split(map_df[[sample_id_col]], map_df[[group_col]])
  group_cols <- group_cols[sapply(group_cols, length) > 0]
  if (length(group_cols) == 0) {
    stop("没有任何可用分组（请检查 sample_id_col 是否能匹配 mat 的列名）。")
  }

  # =========================
  # 2) 计算每组 CV（按行/蛋白）
  # =========================
  cv_list <- lapply(names(group_cols), function(g) {
    cols_g <- group_cols[[g]]
    subm <- mat[, cols_g, drop = FALSE]
    n_g <- ncol(subm)

    # 每组需要的最少非 NA 数：max(min_non_na, ceil(n_g * detect_rate))
    need_non_na <- max(min_non_na, ceiling(n_g * min_detect_rate))

    apply(subm, 1, function(x) {
      x <- as.numeric(x)
      nn <- sum(!is.na(x))
      if (nn < need_non_na) {
        return(NA_real_)
      }

      m <- if (na_rm) mean(x, na.rm = TRUE) else mean(x)
      s <- if (na_rm) sd(x, na.rm = TRUE) else sd(x)

      if (is.na(m) || is.na(s) || abs(m) < min_mean) {
        return(NA_real_)
      }
      (s / m) * 100
    })
  })
  names(cv_list) <- names(group_cols)

  cv_wide <- as.data.frame(cv_list, check.names = FALSE)

  cv_long <- cv_wide |>
    tibble::rownames_to_column(var = "Protein") |>
    tidyr::pivot_longer(
      cols = -Protein,
      names_to = "Group",
      values_to = "CV"
    ) |>
    tidyr::drop_na(CV)

  # =========================
  # 3) 组顺序（Group factor levels）
  # =========================
  if (!is.null(group_levels)) {
    group_levels_use <- as.character(group_levels)
    cv_long$Group <- factor(cv_long$Group, levels = group_levels_use)
  } else {
    group_levels_use <- names(group_cols)
    cv_long$Group <- factor(cv_long$Group, levels = group_levels_use)
  }

  # =========================
  # 4) CV cap（只用于展示）
  # =========================
  if (!is.null(cv_cap)) {
    cv_long$CV_plot <- pmin(cv_long$CV, cv_cap)
    y_use <- "CV_plot"
  } else {
    y_use <- "CV"
  }

  # =========================
  # 5) ★新增：图例标签加上每组样本量 n
  #    这里的 n 是“每组样本数”（来自 group_cols），不是 CV 行数
  # =========================
  n_vec <- vapply(
    group_levels_use,
    function(g) {
      if (!is.null(group_cols[[g]])) length(group_cols[[g]]) else 0L
    },
    integer(1)
  )

  legend_labels <- stats::setNames(
    paste0(group_levels_use, " (n=", n_vec, ")"),
    group_levels_use
  )

  # =========================
  # 6) 画图
  # =========================
  p <- ggplot2::ggplot(
    cv_long,
    ggplot2::aes(x = Group, y = .data[[y_use]], fill = Group)
  )

  ## 1. 最底层：小提琴图
  if (isTRUE(add_violin)) {
    p <- p +
      ggplot2::geom_violin(
        trim = FALSE,
        alpha = 0.6
      )
  }

  ## 2. 中间层：散点
  if (isTRUE(add_points)) {
    p <- p +
      ggplot2::geom_jitter(
        shape = 21,
        color = if (isTRUE(point_border)) "black" else NA,
        stroke = if (isTRUE(point_border)) point_stroke else 0,
        width = point_width,
        size = point_size,
        alpha = point_alpha
      )
  }

  ## 3. 最上层：箱线图
  if (isTRUE(add_box)) {
    p <- p +
      ggplot2::geom_boxplot(
        width = 0.12,
        outlier.shape = NA,
        alpha = 0.4
      )
  }

  bg_fill <- if (isTRUE(transparent_bg)) "transparent" else "white"

  p <- p +
    ggplot2::labs(
      title = title %||% "",
      subtitle = subtitle,
      y = ylab,
      x = NULL,
      fill = group_col # 图例标题用 group_col
    ) +
    ggplot2::theme_minimal(base_family = font_family) +
    ggplot2::theme(
      # ★原来是 none，这里打开图例（你要显示 n）
      legend.position = "right",

      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),

      plot.title = ggplot2::element_text(
        size = fs_title,
        face = if (isTRUE(bold_title)) "bold" else "plain",
        hjust = 0.5
      ),
      plot.subtitle = ggplot2::element_text(size = fs_subtitle, hjust = 0.5),

      axis.title.x = ggplot2::element_text(
        size = fs_axis_title,
        face = if (isTRUE(bold_axis_title)) "bold" else "plain"
      ),
      axis.title.y = ggplot2::element_text(
        size = fs_axis_title,
        face = if (isTRUE(bold_axis_title)) "bold" else "plain"
      ),
      axis.text.x = ggplot2::element_text(
        size = fs_axis_text,
        angle = angle,
        hjust = hjust,
        vjust = vjust,
        color = "black"
      ),
      axis.text.y = ggplot2::element_text(size = fs_axis_text, color = "black"),

      strip.text = ggplot2::element_text(
        size = fs_strip_text,
        face = if (isTRUE(bold_strip_text)) "bold" else "plain"
      ),

      panel.border = ggplot2::element_rect(
        color = "black",
        fill = NA,
        linewidth = panel_border_size
      ),

      panel.background = ggplot2::element_rect(fill = bg_fill, colour = NA),
      plot.background = ggplot2::element_rect(fill = bg_fill, colour = NA),
      legend.background = ggplot2::element_rect(fill = bg_fill, colour = NA),
      legend.box.background = ggplot2::element_rect(fill = bg_fill, colour = NA)
    )

  # =========================
  # 7) 顶部 Median 标签
  # =========================
  top_df <- cv_long |>
    dplyr::group_by(Group) |>
    dplyr::summarise(
      med = median(CV, na.rm = TRUE),
      y_top = max(.data[[y_use]], na.rm = TRUE),
      .groups = "drop"
    )

  y_range <- diff(range(cv_long[[y_use]], na.rm = TRUE))
  if (is.na(y_range) || y_range == 0) {
    y_range <- 1
  }
  top_df <- top_df |>
    dplyr::mutate(y_top = y_top + 0.3 * y_range)

  p <- p +
    ggplot2::geom_text(
      data = top_df,
      ggplot2::aes(
        x = Group,
        y = Inf - 0.1,
        label = paste0("Median=", sprintf("%.2f", med))
      ),
      inherit.aes = FALSE,
      vjust = 2,
      size = 4,
      fontface = if (isTRUE(bold_median_label)) "bold" else "plain",
      color = "black"
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0.02, 0.12))
    )

  # =========================
  # 8) 颜色 + ★图例 labels（带 n）
  # =========================
  if (!is.null(cols)) {
    cols <- as.character(cols)
    has_names <- !is.null(names(cols)) && any(nzchar(names(cols)))

    if (has_names) {
      miss <- setdiff(group_levels_use, names(cols))
      if (length(miss) > 0) {
        stop(
          "cols 缺少这些组的颜色：",
          paste(miss, collapse = ", "),
          "\n当前 group_levels: ",
          paste(group_levels_use, collapse = ", ")
        )
      }
      cols_use <- cols[group_levels_use]
    } else {
      if (length(cols) < length(group_levels_use)) {
        stop(
          "cols 颜色数量不足：给了 ",
          length(cols),
          " 个，但需要 ",
          length(group_levels_use),
          " 个。"
        )
      }
      cols_use <- cols[seq_along(group_levels_use)]
      names(cols_use) <- group_levels_use
    }

    p <- p +
      ggplot2::scale_fill_manual(
        values = cols_use,
        breaks = group_levels_use,
        labels = legend_labels[group_levels_use],
        drop = FALSE
      )
  } else {
    p <- p +
      ggplot2::scale_fill_discrete(
        breaks = group_levels_use,
        labels = legend_labels[group_levels_use],
        drop = FALSE
      )
  }

  # =========================
  # 9) cap 提示（可选）
  # =========================
  if (!is.null(cv_cap) && is.null(subtitle)) {
    p <- p +
      ggplot2::labs(
        subtitle = paste0("CV capped at ", cv_cap, " for visualization")
      )
  }

  # =========================
  # 10) 保存（可选）
  # =========================
  if (!is.null(out_prefix)) {
    out_dir <- dirname(out_prefix)
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }

    save_formats <- tolower(save_formats)
    bg_arg <- if (isTRUE(transparent_bg)) "transparent" else "white"

    if ("pdf" %in% save_formats) {
      ggplot2::ggsave(
        paste0(out_prefix, ".pdf"),
        plot = p,
        width = width,
        height = height,
        bg = bg_arg
      )
    }
    if ("png" %in% save_formats) {
      ggplot2::ggsave(
        paste0(out_prefix, ".png"),
        plot = p,
        width = width,
        height = height,
        dpi = dpi,
        bg = bg_arg
      )
    }
    if ("tiff" %in% save_formats || "tif" %in% save_formats) {
      ggplot2::ggsave(
        paste0(out_prefix, ".tiff"),
        plot = p,
        width = width,
        height = height,
        dpi = dpi,
        bg = bg_arg,
        compression = tiff_compression
      )
    }
  }

  return(list(
    cv_long = cv_long,
    cv_wide = cv_wide,
    plot = p,
    group_cols = group_cols,
    legend_labels = legend_labels # ★返回一下，方便你检查各组 n
  ))
}

##########################################################
###################### 新增 鉴定量 CV #####################
library(tidyverse)

plot_identification_count_cv_violin <- function(
  mat,
  metadata,
  sample_id_col = "sample_id",
  group_col = "group",
  group_levels = NULL,
  cols = NULL,
  zero_as_na = FALSE,

  add_violin = TRUE,
  add_box = TRUE,
  add_points = FALSE,

  point_alpha = 0.8,
  point_size = 1.8,
  point_width = 0.15,
  point_border = TRUE,
  point_stroke = 0.1,

  title = NULL,
  subtitle = NULL,
  ylab = "Identification Count",
  font_family = NULL,
  fs_title = 16,
  fs_subtitle = 12,
  fs_axis_title = 14,
  fs_axis_text = 12,
  panel_border_size = 1,
  transparent_bg = FALSE,

  bold_title = FALSE,
  bold_axis_title = FALSE,
  bold_cv_label = FALSE,

  angle = NULL,
  hjust = 0.5,
  vjust = NULL,

  out_prefix = NULL,
  save_formats = NULL,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300,
  tiff_compression = "lzw"
) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  ## =========================
  ## 0) 输入检查
  ## =========================
  if (is.null(colnames(mat))) {
    stop("mat 没有列名 colnames(mat)，列名必须是样本名。")
  }

  mat <- as.matrix(mat)

  if (!is.data.frame(metadata)) {
    stop("metadata 必须是 data.frame。")
  }

  if (!sample_id_col %in% colnames(metadata)) {
    stop(sprintf("metadata 中找不到 sample_id_col='%s'。", sample_id_col))
  }

  if (!group_col %in% colnames(metadata)) {
    stop(sprintf("metadata 中找不到 group_col='%s'。", group_col))
  }

  if (isTRUE(zero_as_na)) {
    mat[mat == 0] <- NA
  }

  ## =========================
  ## 1) 建立 sample -> group 映射
  ## =========================
  map_df <- metadata %>%
    transmute(
      Sample = as.character(.data[[sample_id_col]]),
      Group = as.character(.data[[group_col]])
    )

  map_df <- map_df %>%
    filter(Sample %in% colnames(mat))

  miss_samples <- setdiff(colnames(mat), map_df$Sample)

  if (length(miss_samples) > 0) {
    warning(
      "这些样本在 metadata 中找不到，将不会参与鉴定数 CV 计算：\n",
      paste(miss_samples, collapse = ", ")
    )
  }

  if (nrow(map_df) == 0) {
    stop("mat 的列名和 metadata 的样本名没有匹配上。")
  }

  ## =========================
  ## 2) 统计每个样本的鉴定量
  ## =========================
  ident_df <- tibble(
    Sample = colnames(mat),
    Ident_Count = colSums(!is.na(mat))
  ) %>%
    left_join(map_df, by = "Sample") %>%
    filter(!is.na(Group), Group != "")

  ## =========================
  ## 3) 组顺序
  ## =========================
  if (!is.null(group_levels)) {
    group_levels_use <- as.character(group_levels)
  } else {
    group_levels_use <- unique(ident_df$Group)
  }

  ident_df$Group <- factor(ident_df$Group, levels = group_levels_use)

  ## =========================
  ## 4) 每组计算鉴定量 CV
  ## =========================
  stat_df <- ident_df %>%
    group_by(Group) %>%
    summarise(
      n = sum(!is.na(Ident_Count)),
      mean_ident = mean(Ident_Count, na.rm = TRUE),
      sd_ident = sd(Ident_Count, na.rm = TRUE),
      cv_percent = ifelse(
        n >= 2 & !is.na(mean_ident) & mean_ident != 0,
        sd_ident / mean_ident * 100,
        NA_real_
      ),
      .groups = "drop"
    )

  ## 图例标签：组名 + 样本数
  legend_labels <- stats::setNames(
    paste0(stat_df$Group, " (n=", stat_df$n, ")"),
    stat_df$Group
  )

  ## =========================
  ## 5) 颜色
  ## =========================
  if (!is.null(cols)) {
    cols <- as.character(cols)
    has_names <- !is.null(names(cols)) && any(nzchar(names(cols)))

    if (has_names) {
      miss <- setdiff(group_levels_use, names(cols))

      if (length(miss) > 0) {
        stop(
          "cols 缺少这些组的颜色：",
          paste(miss, collapse = ", ")
        )
      }

      cols_use <- cols[group_levels_use]
    } else {
      if (length(cols) < length(group_levels_use)) {
        stop(
          "cols 颜色数量不足：给了 ",
          length(cols),
          " 个，但需要 ",
          length(group_levels_use),
          " 个。"
        )
      }

      cols_use <- cols[seq_along(group_levels_use)]
      names(cols_use) <- group_levels_use
    }
  } else {
    cols_use <- setNames(
      scales::hue_pal()(length(group_levels_use)),
      group_levels_use
    )
  }

  ## =========================
  ## 6) 画小提琴图
  ## =========================
  p <- ggplot(
    ident_df,
    aes(
      x = Group,
      y = Ident_Count,
      fill = Group
    )
  )

  if (isTRUE(add_violin)) {
    p <- p +
      geom_violin(trim = FALSE, alpha = 0.6)
  }

  if (isTRUE(add_box)) {
    p <- p +
      geom_boxplot(
        width = 0.12,
        outlier.shape = NA,
        alpha = 0.4
      )
  }

  if (isTRUE(add_points)) {
    p <- p +
      geom_jitter(
        shape = 21,
        color = if (isTRUE(point_border)) "black" else NA,
        stroke = if (isTRUE(point_border)) point_stroke else 0,
        width = point_width,
        size = point_size,
        alpha = point_alpha
      )
  }

  bg_fill <- if (isTRUE(transparent_bg)) "transparent" else "white"

  p <- p +
    labs(
      title = title %||% "",
      subtitle = subtitle,
      y = ylab,
      x = NULL,
      fill = group_col
    ) +
    theme_minimal(base_family = font_family) +
    theme(
      legend.position = "right",

      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),

      plot.title = element_text(
        size = fs_title,
        face = if (isTRUE(bold_title)) "bold" else "plain",
        hjust = 0.5
      ),
      plot.subtitle = element_text(
        size = fs_subtitle,
        hjust = 0.5
      ),

      axis.title.x = element_text(
        size = fs_axis_title,
        face = if (isTRUE(bold_axis_title)) "bold" else "plain"
      ),
      axis.title.y = element_text(
        size = fs_axis_title,
        face = if (isTRUE(bold_axis_title)) "bold" else "plain"
      ),

      axis.text.x = element_text(
        size = fs_axis_text,
        angle = angle,
        hjust = hjust,
        vjust = vjust,
        color = "black"
      ),
      axis.text.y = element_text(
        size = fs_axis_text,
        color = "black"
      ),

      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = panel_border_size
      ),

      panel.background = element_rect(fill = bg_fill, colour = NA),
      plot.background = element_rect(fill = bg_fill, colour = NA),
      legend.background = element_rect(fill = bg_fill, colour = NA),
      legend.box.background = element_rect(fill = bg_fill, colour = NA)
    )

  ## =========================
  ## 7) 顶部 CV 标签
  ## 关键：和强度 CV 图一致，固定画在最上方
  ## =========================
  p <- p +
    geom_text(
      data = stat_df,
      aes(
        x = Group,
        y = Inf - 0.1,
        label = paste0("CV=", sprintf("%.2f", cv_percent), "%")
      ),
      inherit.aes = FALSE,
      vjust = 2,
      size = 4,
      fontface = if (isTRUE(bold_cv_label)) "bold" else "plain",
      color = "black"
    ) +
    coord_cartesian(clip = "off") +
    scale_y_continuous(
      expand = expansion(mult = c(0.02, 0.12))
    )

  ## =========================
  ## 8) 颜色和图例
  ## =========================
  p <- p +
    scale_fill_manual(
      values = cols_use,
      breaks = group_levels_use,
      labels = legend_labels[group_levels_use],
      drop = FALSE
    )

  ## =========================
  ## 9) 保存
  ## =========================
  if (!is.null(out_prefix)) {
    out_dir <- dirname(out_prefix)

    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }

    save_formats <- tolower(save_formats)
    bg_arg <- if (isTRUE(transparent_bg)) "transparent" else "white"

    if ("pdf" %in% save_formats) {
      ggsave(
        paste0(out_prefix, ".pdf"),
        plot = p,
        width = width,
        height = height,
        units = units,
        bg = bg_arg
      )
    }

    if ("png" %in% save_formats) {
      ggsave(
        paste0(out_prefix, ".png"),
        plot = p,
        width = width,
        height = height,
        units = units,
        dpi = dpi,
        bg = bg_arg
      )
    }

    if ("tiff" %in% save_formats || "tif" %in% save_formats) {
      ggsave(
        paste0(out_prefix, ".tiff"),
        plot = p,
        width = width,
        height = height,
        units = units,
        dpi = dpi,
        bg = bg_arg,
        compression = tiff_compression
      )
    }
  }

  return(list(
    sample_ident_df = ident_df,
    group_stat_df = stat_df,
    plot = p,
    legend_labels = legend_labels
  ))
}

if (F) {
  res_protein_ident <- plot_identification_count_cv_violin(
    mat = matrix_all,
    metadata = metadata2,
    sample_id_col = "sample_id",
    group_col = "group",
    title = "Protein Identification Count",
    ylab = "Protein Identification Count",
    add_points = TRUE,
    cols = c(
      "251213-GO-007" = "#4E79A7",
      "260116-GO-007" = "#F28E2B",
      "260206-GO-007" = "#59A14F",
      "260210-GO-007" = "#E15759"
    ),
    out_prefix = "result/protein_identification_count",
    save_formats = c("pdf", "png")
  )

  res_protein_ident$plot
}
