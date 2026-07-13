plot_peptide_cv_box <- function(mat_expr,
                                value_col = "ratio_LH",   # ★ 数值列名（pivot_longer 输出的 values_to）
                                peptide_col = "肽段",
                                rep_col    = "重复测定",
                                min_non_na = 3,
                                min_mean   = 1e-8,
                                jitter_width = 0.15,
                                jitter_alpha = 0.35,
                                jitter_size  = 0.8,
                                box_width    = 0.6,
                                xlab = "肽段（按CV从小到大排序）",
                                ylab = "Light intensity",
                                title = "每个肽段的箱线图 + CV（仅标注最小/中间/最大）",
                                rotate_x_text = TRUE,
                                ylim_expand_frac = 0.12,
                                label_inset_frac = 0.03,
                                out_file = NULL,# ★建议传“完整文件名”，如 "xxx.png"
                                width = 10,
                                height = 6,
                                dpi = 300,
                                max_cv = 20) {
  
  # ----------------------------
  # 依赖包检查
  # ----------------------------
  stopifnot(!missing(mat_expr))
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("需要安装 dplyr")
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("需要安装 tidyr")
  if (!requireNamespace("tibble", quietly = TRUE)) stop("需要安装 tibble")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("需要安装 ggplot2")
  
  # 非标准求值用
  val_sym <- rlang::sym(value_col)
  pep_sym <- rlang::sym(peptide_col)
  rep_sym <- rlang::sym(rep_col)
  
  # ----------------------------
  # 1) 宽表 -> 长表
  # ----------------------------
  df_wide <- as.data.frame(mat_expr, check.names = FALSE)
  if (is.null(rownames(df_wide))) stop("mat_expr 必须有行名（肽段）作为 rownames。")
  
  sample_long <- df_wide %>%
    tibble::rownames_to_column("Peptides") %>%
    tidyr::pivot_longer(
      cols = -Peptides,
      names_to  = rep_col,
      values_to = value_col
    ) %>%
    dplyr::transmute(
      !!peptide_col := as.character(.data$Peptides),
      !!rep_col     := as.character(.data[[rep_col]]),
      !!value_col   := suppressWarnings(as.numeric(.data[[value_col]]))
    ) %>%
    dplyr::filter(!is.na(!!val_sym))
  
  # ----------------------------
  # 2) 计算每个肽段 CV
  # ----------------------------
  cv_df <- sample_long %>%
    dplyr::group_by(!!pep_sym) %>%
    dplyr::summarise(
      n        = sum(!is.na(!!val_sym)),
      mean_val = mean(!!val_sym, na.rm = TRUE),
      sd_val   = stats::sd(!!val_sym, na.rm = TRUE),
      CV = dplyr::if_else(
        n < min_non_na | is.na(mean_val) | abs(mean_val) < min_mean,
        NA_real_,
        sd_val / mean_val * 100
      ),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      label = dplyr::if_else(is.na(.data$CV), "CV=NA", sprintf("CV=%.1f%%", .data$CV))
    )
  
  # ----------------------------
  # 3) CV 回写到长表
  # ----------------------------
  df_long2 <- sample_long %>%
    dplyr::left_join(
      cv_df %>% dplyr::select(!!pep_sym, .data$CV, .data$label),
      by = rlang::set_names(peptide_col, peptide_col)
    )
  
  # ----------------------------
  # 4) 按 CV 排序（NA 放最后），用于 x 轴顺序
  # ----------------------------
  peptide_order <- cv_df %>%
    dplyr::arrange(is.na(.data$CV), .data$CV) %>%
    dplyr::pull(!!pep_sym)
  
  df_long2 <- df_long2 %>%
    dplyr::mutate(!!pep_sym := factor(!!pep_sym, levels = peptide_order))
  
  cv_df <- cv_df %>%
    dplyr::mutate(!!pep_sym := factor(!!pep_sym, levels = peptide_order))
  
  # ----------------------------
  # 5) 只取最小/最大/中间 CV 三个用于标注
  # ----------------------------
  cv_non_na <- cv_df %>%
    dplyr::filter(!is.na(.data$CV)) %>%
    dplyr::arrange(.data$CV)
  
  n_cv <- nrow(cv_non_na)
  if (n_cv == 0) {
    cv_anno <- cv_non_na
  } else {
    mid_idx <- ceiling(n_cv / 2)
    cv_anno <- dplyr::bind_rows(
      cv_non_na %>% dplyr::slice(1),
      cv_non_na %>% dplyr::slice(n_cv),
      cv_non_na %>% dplyr::slice(mid_idx)
    ) %>%
      dplyr::distinct(!!pep_sym, .keep_all = TRUE)
  }
  
  # ----------------------------
  # 6) 统一标签 y 位置
  # ----------------------------
  y_max <- max(df_long2[[value_col]], na.rm = TRUE)
  y_min <- min(df_long2[[value_col]], na.rm = TRUE)
  y_rng <- y_max - y_min
  if (!is.finite(y_rng) || y_rng == 0) y_rng <- 1
  
  y_upper_plot <- y_max + ylim_expand_frac * y_rng
  y_inset      <- label_inset_frac * y_rng
  
  if (nrow(cv_anno) > 0) {
    cv_anno <- cv_anno %>%
      dplyr::mutate(y_pos = y_upper_plot - y_inset)
  }
  
  # ----------------------------
  # 7) 作图
  # ----------------------------
  p <- ggplot2::ggplot(df_long2, ggplot2::aes(x = !!pep_sym, y = !!val_sym)) +
    ggplot2::geom_boxplot(width = box_width, outlier.shape = NA) +
    ggplot2::geom_jitter(width = jitter_width, alpha = jitter_alpha, size = jitter_size)
  
  if (nrow(cv_anno) > 0) {
    p <- p +
      ggplot2::geom_text(
        data = cv_anno,
        ggplot2::aes(x = !!pep_sym, y = .data$y_pos, label = .data$label),
        inherit.aes = FALSE,
        size = 3,
        vjust = 0,
        fontface = "bold"
      )
  }
  
  p <- p +
    ggplot2::coord_cartesian(ylim = c(NA, y_upper_plot), clip = "off") +
    ggplot2::labs(x = xlab, y = ylab, title = title) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title  = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text.x = if (rotate_x_text) ggplot2::element_text(angle = 45, hjust = 1) else ggplot2::element_text()
    )
  # 保存图片（如果 out_file 不为空）
  if (!is.null(out_file) && nzchar(out_file)) {
    dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(filename = out_file, plot = p, width = width, height = height, dpi = dpi)
  }
  # ----------------------------
  # 8) 只保留“肽段-CV”对照表
  # ----------------------------
  peptide_cv_tbl <- cv_df %>%
    dplyr::transmute(
      !!peptide_col := as.character(!!pep_sym),
      CV = as.numeric(.data$CV)
    ) %>%
    dplyr::arrange(is.na(.data$CV), .data$CV)
  
  return(list(
    plot = p,
    df_long = df_long2,
    cv_df = cv_df,
    cv_anno = cv_anno,
    peptide_cv_tbl = as.data.frame(peptide_cv_tbl),
    filtered_peptide = as.data.frame(peptide_cv_tbl[peptide_cv_tbl$CV <= max_cv,1],drop = F)
  ))
}

if(F) {
  res <- plot_peptide_cv_box(light_pics_expr, value_col = "intensity")
  
}