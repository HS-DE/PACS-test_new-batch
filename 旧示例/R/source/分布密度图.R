plot_plate_density <- function(mat,
                               metadata,
                               sample_id_col = "sample_id",
                               plate_col = "Plate",
                               color_by = c("plate", "sample", "sample_by_plate"),
                               cols = NULL,
                               color_alpha = 1,          # ★颜色本身透明度（0~1）
                               legend = TRUE,
                               log2_transform = TRUE,
                               add1 = 1,
                               mode = c("values", "sample_median"),
                               n_per_plate = NULL,
                               title = NULL,
                               subtitle = NULL,
                               font_family = NULL,
                               fs_title = 16,
                               fs_subtitle = 12,
                               fs_axis_title = 14,
                               fs_axis_text = 12,
                               fs_strip_text = 20,
                               bold_title = FALSE,
                               bold_axis_title = FALSE,
                               bold_strip = FALSE,
                               panel_border_size = 1,
                               linewidth = 0.5,
                               caption = NULL,
                               transparent_bg = FALSE) {
  
  mode <- match.arg(mode)
  color_by <- match.arg(color_by)
  stopifnot(!is.null(colnames(mat)))
  mat <- as.matrix(mat)
  
  # color_alpha 检查
  color_alpha <- as.numeric(color_alpha)
  if (!is.finite(color_alpha) || color_alpha < 0 || color_alpha > 1) {
    stop("color_alpha 必须是 0~1 之间的数值。")
  }
  
  if (!is.data.frame(metadata)) stop("metadata 必须是 data.frame。")
  if (!sample_id_col %in% colnames(metadata)) stop(sprintf("metadata 中找不到 sample_id_col='%s'。", sample_id_col))
  if (!plate_col %in% colnames(metadata)) stop(sprintf("metadata 中找不到 plate_col='%s'。", plate_col))
  
  # ★保留用户指定的 Plate 顺序
  plate_levels0 <- if (is.factor(metadata[[plate_col]])) levels(metadata[[plate_col]]) else NULL
  
  meta2 <- metadata[metadata[[sample_id_col]] %in% colnames(mat), , drop = FALSE]
  meta2[[sample_id_col]] <- as.character(meta2[[sample_id_col]])
  # ★不要直接 as.character，保留 levels
  if (!is.null(plate_levels0)) {
    meta2[[plate_col]] <- factor(as.character(meta2[[plate_col]]), levels = plate_levels0)
  } else {
    meta2[[plate_col]] <- as.character(meta2[[plate_col]])
  }
  
  
  # =========================
  # 1) 组装 df
  # =========================
  if (mode == "sample_median") {
    x <- apply(mat[, meta2[[sample_id_col]], drop = FALSE], 2, median, na.rm = TRUE)
    df <- data.frame(sample = names(x), value = as.numeric(x), stringsAsFactors = FALSE)
    df <- dplyr::left_join(df, meta2, by = setNames(sample_id_col, "sample"))
  } else {
    df <- data.frame(mat, check.names = FALSE) |>
      tibble::rownames_to_column("Feature") |>
      tidyr::pivot_longer(cols = -Feature, names_to = "sample", values_to = "value") |>
      tidyr::drop_na(value)
    
    df <- dplyr::left_join(df, meta2, by = setNames(sample_id_col, "sample"))
    # ★保证 df 里的 Plate 也按用户指定 levels
    if (!is.null(plate_levels0)) {
      df[[plate_col]] <- factor(as.character(df[[plate_col]]), levels = plate_levels0)
    }
    if (!is.null(n_per_plate)) {
      df <- df %>%
        dplyr::group_by(.data[[plate_col]]) %>%
        dplyr::group_modify(~{
          n_take <- min(n_per_plate, nrow(.x))
          dplyr::slice_sample(.x, n = n_take)
        }) %>%
        dplyr::ungroup()
    }
  }
  
  if (log2_transform) df$value <- log2(df$value + add1)
  if (is.null(title)) title <- paste0("Density (mode=", mode, ")")
  xlab <- if (log2_transform) paste0("log2(Intensity+", add1, ")") else "Intensity"
  
  # =========================
  # 2) 映射颜色变量
  # =========================
  if (color_by == "sample") {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = value, group = sample, colour = sample))
    color_var <- "sample"
    legend_title <- "Sample"
  } else if (color_by == "sample_by_plate") {
    p <- ggplot2::ggplot(df, ggplot2::aes(
      x = value, group = sample, colour = .data[[plate_col]]
    ))
    color_var <- plate_col
    legend_title <- plate_col
  } else { # "plate"
    p <- ggplot2::ggplot(df, ggplot2::aes(
      x = value, colour = .data[[plate_col]]
    ))
    color_var <- plate_col
    legend_title <- plate_col
  }
  
  # =========================
  # 3) ★图例标签：每组追加 n=样本量
  #    - plate: n = 该 plate 下 distinct(sample) 数
  #    - sample: n = 1
  # =========================
  if (identical(color_var, "sample")) {
    n_tbl <- df %>%
      dplyr::distinct(sample) %>%
      dplyr::count(sample, name = "n")  # 每个 sample 基本都是 1
    n_tbl$lv <- as.character(n_tbl$sample)
  } else {
    n_tbl <- df %>%
      dplyr::distinct(sample, .data[[plate_col]]) %>%
      dplyr::count(.data[[plate_col]], name = "n")
    n_tbl$lv <- as.character(n_tbl[[plate_col]])
  }
  
  # 当前图例水平（按出现的水平排序，保证 labels 对得上）
  if (identical(color_var, "sample")) {
    levs <- unique(as.character(df$sample))   # sample 这里你想按列顺序的话也可以改成 colnames(mat)
  } else {
    # ★优先用用户指定的 Plate levels
    if (!is.null(plate_levels0)) {
      levs <- plate_levels0
    } else if (is.factor(df[[plate_col]])) {
      levs <- levels(df[[plate_col]])
    } else {
      levs <- unique(as.character(df[[plate_col]]))
    }
  }
  levs <- levs[!is.na(levs)]
  
  
  n_map <- stats::setNames(n_tbl$n, n_tbl$lv)
  n_vec <- n_map[levs]
  n_vec[is.na(n_vec)] <- 0L
  
  legend_labels <- stats::setNames(
    paste0(levs, " (n=", n_vec, ")"),
    levs
  )
  
  # =========================
  # 4) 主题
  # =========================
  bg_fill <- if (isTRUE(transparent_bg)) "transparent" else "white"
  
  p <- p +
    ggplot2::geom_density(linewidth = linewidth) +
    ggplot2::labs(
      title = title, subtitle = subtitle, caption = caption,
      x = xlab, y = "Density",
      colour = legend_title
    ) +
    ggplot2::theme_bw(base_family = font_family) +
    ggplot2::theme(
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
      
      axis.text.x = ggplot2::element_text(size = fs_axis_text, hjust = 0.5, color = "black"),
      axis.text.y = ggplot2::element_text(size = fs_axis_text, color = "black"),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = panel_border_size),
      
      panel.background = ggplot2::element_rect(fill = bg_fill, colour = NA),
      plot.background  = ggplot2::element_rect(fill = bg_fill, colour = NA),
      legend.background = ggplot2::element_rect(fill = bg_fill, colour = NA),
      legend.box.background = ggplot2::element_rect(fill = bg_fill, colour = NA)
    )
  
  # =========================
  # 5) 只改“颜色本身透明度” + 颜色尺度 + ★labels(带 n)
  # =========================
  add_color_alpha <- function(col_vec, a) {
    if (isTRUE(all.equal(a, 1))) return(col_vec)
    if (requireNamespace("scales", quietly = TRUE)) {
      scales::alpha(col_vec, a)
    } else {
      grDevices::adjustcolor(col_vec, alpha.f = a)
    }
  }
  
  # 需要 manual scale 的条件：
  # - 用户给了 cols
  # - 或者 color_alpha < 1（需要把 alpha 写进颜色）
  need_manual_scale <- !is.null(cols) || !isTRUE(all.equal(color_alpha, 1))
  
  if (need_manual_scale) {
    
    if (is.null(cols)) {
      # 自动生成颜色
      if (requireNamespace("scales", quietly = TRUE)) {
        cols_use <- scales::hue_pal()(length(levs))
      } else {
        cols_use <- grDevices::hcl.colors(length(levs), palette = "Dynamic")
      }
      names(cols_use) <- levs
    } else {
      cols <- as.character(cols)
      has_names <- !is.null(names(cols)) && any(nzchar(names(cols)))
      
      if (has_names) {
        miss <- setdiff(levs, names(cols))
        if (length(miss) > 0) stop("cols 缺少这些水平的颜色：", paste(miss, collapse = ", "))
        cols_use <- cols[levs]
      } else {
        if (length(cols) < length(levs)) {
          stop("cols 颜色数量不足：给了 ", length(cols), " 个，但需要 ", length(levs), " 个。")
        }
        cols_use <- cols[seq_along(levs)]
        names(cols_use) <- levs
      }
    }
    
    cols_use <- add_color_alpha(cols_use, color_alpha)
    
    p <- p + ggplot2::scale_colour_manual(
      values = cols_use,
      breaks = levs,
      labels = legend_labels[levs],
      drop = FALSE
    )
    
  } else {
    # 不手动给颜色，但仍然要改 legend labels（带 n）
    p <- p + ggplot2::scale_colour_discrete(
      breaks = levs,
      labels = legend_labels[levs],
      drop = FALSE
    )
  }
  
  if (!isTRUE(legend)) {
    p <- p + ggplot2::theme(legend.position = "none")
  }
  
  p
}
