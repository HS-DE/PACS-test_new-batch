plot_cv_multi_stage <- function(
    mat_list,
    metadata,
    sample_id_col = "sample_id",
    group_col = "test_group",
    group_levels = NULL,
    stage_levels = NULL,
    cols_stage = NULL,
    min_detect_rate = 0.70,
    min_non_na = 3,
    cv_cap = NULL,
    title = NULL,
    ylab = "Coefficient of Variation (%)",
    fs_axis_text = 12,
    fs_axis_title = 14,
    fs_title = 16,
    bold_title = FALSE,
    bold_axis_title = FALSE,
    bold_median_label = FALSE,
    
    add_violin = TRUE,
    add_box = TRUE,
    dodge_width = 0.8,
    label_median = TRUE,
    median_fmt = "%.2f",
    median_prefix = "Median=",
    median_y_mode = c("per_group_top", "panel_top"),
    median_y_pad = 0.06,
    alpha_violin = 0.7,
    alpha_box = 0.7,
    fs_legend_title = 14,
    fs_legend_text  = 14,
    
    # ★只在横坐标显示每组 n
    show_group_n_on_x = TRUE,
    
    out_file = NULL,
    width = 8,
    height = 5,
    dpi = 300
) {
  median_y_mode <- match.arg(median_y_mode)
  source("C:/Work/SH/code/source/CV作图.R")
  # ---------- 0) 检查 ----------
  if (!is.list(mat_list) || length(mat_list) < 2) {
    stop("mat_list 必须是一个长度>=2的 list，并且建议是命名 list：list(Stage1=mat1, Stage2=mat2, ...)")
  }
  if (is.null(names(mat_list)) || any(names(mat_list) == "")) {
    stop("mat_list 必须是命名 list，每个元素名将作为 Stage 名。")
  }
  if (!is.data.frame(metadata)) stop("metadata 必须是 data.frame。")
  if (!sample_id_col %in% colnames(metadata)) stop(sprintf("metadata 中找不到 sample_id_col='%s'。", sample_id_col))
  if (!group_col %in% colnames(metadata)) stop(sprintf("metadata 中找不到 group_col='%s'。", group_col))
  
  # stage 顺序
  if (is.null(stage_levels)) stage_levels <- names(mat_list)
  stage_levels <- as.character(stage_levels)
  if (!all(stage_levels %in% names(mat_list))) {
    stop("stage_levels 必须都在 names(mat_list) 里。")
  }
  
  # group 顺序
  if (is.null(group_levels)) {
    group_levels <- levels(factor(metadata[[group_col]]))
  }
  group_levels <- as.character(group_levels)
  
  # ---------- 1) 计算每个 Group 的样本量 n（用于横坐标标签） ----------
  meta_map <- metadata[, c(sample_id_col, group_col), drop = FALSE]
  meta_map[[sample_id_col]] <- as.character(meta_map[[sample_id_col]])
  meta_map[[group_col]]     <- as.character(meta_map[[group_col]])
  
  stage_samples <- lapply(stage_levels, function(stg) colnames(mat_list[[stg]]))
  names(stage_samples) <- stage_levels
  union_samples  <- unique(unlist(stage_samples, use.names = FALSE))
  common_samples <- Reduce(intersect, stage_samples)
  
  # 如果不同 stage 样本集合不一致，提醒一下（n 按 union 统计）
  if (length(union_samples) != length(common_samples)) {
    warning("不同 stage 的样本集合不一致：\n",
            "  union_samples = ", length(union_samples), ", common_samples = ", length(common_samples), "\n",
            "横坐标每组 n 按 union_samples（任一 stage 中出现的样本）计算。")
  }
  
  meta_map_use <- meta_map[meta_map[[sample_id_col]] %in% union_samples, , drop = FALSE]
  
  n_df <- meta_map_use |>
    dplyr::distinct(.data[[sample_id_col]], .data[[group_col]]) |>
    dplyr::count(.data[[group_col]], name = "n") |>
    dplyr::mutate(.grp = as.character(.data[[group_col]]))
  
  n_vec <- stats::setNames(n_df$n, n_df$.grp)
  n_vec <- n_vec[group_levels]
  n_vec[is.na(n_vec)] <- 0L
  
  x_labels <- stats::setNames(
    paste0(group_levels, "\n(n=", n_vec, ")"),#paste0("Stage", "\n(n=", n_vec, ")"),#
    group_levels
  )
  
  # ---------- 2) Stage 颜色 ----------
  if (is.null(cols_stage)) {
    cols_stage <- setNames(scales::hue_pal()(length(stage_levels)), stage_levels)
  } else {
    if (is.null(names(cols_stage))) stop("cols_stage 必须是命名向量，names=Stage 名")
    miss <- setdiff(stage_levels, names(cols_stage))
    if (length(miss) > 0) stop("cols_stage 缺少这些 stage 的颜色：", paste(miss, collapse = ", "))
    cols_stage <- cols_stage[stage_levels]
  }
  
  # ---------- 3) 对每个 stage 计算 cv_long ----------
  cv_list <- lapply(stage_levels, function(stg) {
    res <- plot_cv_by_metadata(
      mat = mat_list[[stg]],
      metadata = metadata,
      sample_id_col = sample_id_col,
      group_col = group_col,
      group_levels = group_levels,
      cols = NULL,                   # 组颜色不在这里用（fill 给 Stage）
      min_detect_rate = min_detect_rate,
      min_non_na = min_non_na,
      cv_cap = cv_cap,
      title = NULL
    )
    res$cv_long |>
      dplyr::mutate(Stage = stg)
  })
  cv_all <- dplyr::bind_rows(cv_list)
  
  cv_all$Group <- factor(cv_all$Group, levels = group_levels)
  cv_all$Stage <- factor(cv_all$Stage, levels = stage_levels)
  
  y_use <- if (!is.null(cv_cap) && "CV_plot" %in% colnames(cv_all)) "CV_plot" else "CV"
  pd <- ggplot2::position_dodge(width = dodge_width)
  
  # ---------- 4) 主图 ----------
  p <- ggplot2::ggplot(cv_all, ggplot2::aes(x = Group, y = .data[[y_use]], fill = Stage))
  
  if (isTRUE(add_violin)) {
    p <- p + ggplot2::geom_violin(trim = FALSE, alpha = alpha_violin, position = pd)
  }
  if (isTRUE(add_box)) {
    p <- p + ggplot2::geom_boxplot(width = 0.12, outlier.shape = NA,
                                   alpha = alpha_box, position = pd)
  }
  
  plot_title_el <- if (is.null(title)) {
    ggplot2::element_blank()
  } else {
    ggplot2::element_text(
      size = fs_title,
      face = if (isTRUE(bold_title)) "bold" else "plain",
      hjust = 0.5
    )
  }
  
  p <- p +
    ggplot2::labs(title = title, y = ylab, x = NULL, fill = "Stage") +
    ggplot2::scale_fill_manual(values = cols_stage, drop = FALSE) +
    ggplot2::theme(
      legend.position = "top",
      plot.title = plot_title_el,
      axis.title.y = ggplot2::element_text(
        size = fs_axis_title,
        face = if (isTRUE(bold_axis_title)) "bold" else "plain"
      ),
      axis.text.x  = ggplot2::element_text(size = fs_axis_text, hjust = 0.5, color = "black"),
      axis.text.y  = ggplot2::element_text(size = fs_axis_text, color = "black"),
      legend.title = ggplot2::element_text(size = fs_legend_title),
      legend.text  = ggplot2::element_text(size = fs_legend_text)
    )
  
  # ★只在横坐标显示每组 n
  if (isTRUE(show_group_n_on_x)) {
    p <- p + ggplot2::scale_x_discrete(labels = x_labels, drop = FALSE)
  }
  
  # ---------- 5) 中位数标注（Group + Stage） ----------
  top_df <- NULL
  if (isTRUE(label_median)) {
    top_df <- cv_all |>
      dplyr::group_by(Group, Stage) |>
      dplyr::summarise(
        med = median(CV, na.rm = TRUE),
        y_top = max(.data[[y_use]], na.rm = TRUE),
        .groups = "drop"
      )
    
    if (median_y_mode == "per_group_top") {
      y_range <- diff(range(cv_all[[y_use]], na.rm = TRUE))
      if (is.na(y_range) || y_range == 0) y_range <- 1
      top_df <- top_df |>
        dplyr::mutate(y_lab = y_top + median_y_pad * y_range)
      
      p <- p +
        ggplot2::geom_text(
          data = top_df,
          ggplot2::aes(
            x = Group, y = y_lab,
            label = paste0(median_prefix, sprintf(median_fmt, med)),
            group = Stage
          ),
          position = pd,
          inherit.aes = FALSE,
          vjust = 0,
          size = 4,
          fontface = if (isTRUE(bold_median_label)) "bold" else "plain",
          color = "black"
        ) +
        ggplot2::coord_cartesian(clip = "off") +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.12)))
      
    } else { # panel_top
      p <- p +
        ggplot2::geom_text(
          data = top_df,
          ggplot2::aes(
            x = Group, y = Inf,
            label = paste0(median_prefix, sprintf(median_fmt, med)),
            group = Stage
          ),
          position = pd,
          inherit.aes = FALSE,
          vjust = 1 + median_y_pad,
          size = 4,
          fontface = if (isTRUE(bold_median_label)) "bold" else "plain",
          color = "black"
        ) +
        ggplot2::coord_cartesian(clip = "off") +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.12)))
    }
  }
  
  # ---------- 6) 保存 ----------
  if (!is.null(out_file)) {
    dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(out_file, plot = p, width = width, height = height, dpi = dpi,limitsize = FALSE)
  }
  
  return(list(
    plot = p,
    cv_all = cv_all,
    top_df = top_df,
    y_use = y_use,
    cols_stage = cols_stage,
    stage_levels = stage_levels,
    group_levels = group_levels,
    group_n = n_vec
  ))
}
