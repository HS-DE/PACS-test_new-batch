library(dplyr)
library(tibble)
library(ggplot2)
library(forcats)
library(scales)
plot_detection_rate_rank <- function(
    mat,
    # —— 检出定义：默认非NA算检出 —— #
    detected_fun = function(x) !is.na(x),
    
    # —— x轴设置 —— #
    x_break_by = 100,                 # 每隔多少显示一次rank数字
    x_angle = 90,                     # x轴文字旋转角度（90=逆时针）
    x_text_margin_t = -6,             # x轴文字往图形靠近：负值更近
    x_hjust = 1,
    x_vjust = 0.5,
    
    # —— y轴设置 —— #
    y_limits = c(0, 1),
    y_breaks = 1,                     # 只显示100%
    y_label_accuracy = 1,             # 100% 显示精度
    auto_expand_y_for_labels = TRUE,  # 如果标记线/标签超过1，是否自动放宽y上限
    y_pad_frac = 0.03,                # 顶部留白比例（用于标签不顶到边）
    
    # —— 是否分箱上色 —— #
    color_by_bin = TRUE,
    single_color = "#1b9e77",         # 不分箱时用的统一颜色
    bin_lv = c(">=100%", ">=90%", ">=80%", ">=70%", ">=50%", "<50%"),
    bin_breaks = c(1.00, 0.90, 0.80, 0.70, 0.50),  # 与bin_lv前5项对应
    bin_colors = NULL,                # 不给则用默认配色（见下方）
    
    # —— 图例 / 标题 —— #
    show_legend = TRUE,
    x_lab = "Rank",
    y_lab = "Detection rate",
    
    # —— 注释框 —— #
    show_annotation = TRUE,
    ann_thr = c(1.00, 0.90, 0.80, 0.70, 0.50),
    ann_size = 4,
    ann_label_size = 0.25,
    ann_hjust = 1.02,
    ann_vjust = 1.02,
    
    # —— 外观 —— #
    alpha = 0.9,
    bar_width = 1,
    transparent_bg = TRUE,
    plot_margin = ggplot2::margin(5.5, 80, 5.5, 5.5),  # 右边留白大：用于放注释框
    
    # =========================================================
    # ✅ 标记蛋白（新增功能）
    # =========================================================
    mark_proteins = NULL,             # 需要标记的 Protein 名称向量
    mark_color = "red",               # 标记点颜色
    mark_point_size = 2.2,            # 标记点大小
    mark_line_color = "black",        # 虚线颜色
    mark_line_type = "solid",         # 虚线类型
    mark_text_vjust = 0.5,              # 线到标记基因名的距离
    mark_line_width = 0.8,            # 虚线粗细
    mark_label_size = 4,              # 标签字号（geom_text 的 size）
    mark_label_color = "black",       # 标签颜色
    mark_line_len = 0.05,             # 虚线长度（y尺度单位：det_rate是0~1，所以0.05=5%）
    mark_x_nudge = 0,                 # 所有标签统一x偏移（正=往右）
    mark_x_nudge_map = NULL           # 命名向量：对某些蛋白单独偏移；例如 c(ALB=80)
) {
  # =========================
  # 0) 依赖检查
  # =========================
  pkgs <- c("dplyr", "ggplot2", "scales")
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss) > 0) {
    stop("缺少R包：", paste(miss, collapse = ", "),
         "\n请先安装：install.packages(c('dplyr','ggplot2','scales'))")
  }
  
  # =========================
  # 1) 输入 mat 处理：matrix / data.frame 都可
  # =========================
  if (is.data.frame(mat)) mat <- as.matrix(mat)
  if (!is.matrix(mat)) stop("mat 必须是 matrix 或 data.frame（行=Protein，列=样本）")
  
  prot <- rownames(mat)
  if (is.null(prot)) prot <- paste0("Protein_", seq_len(nrow(mat)))
  
  # =========================
  # 2) 计算检出率 det_rate（按行）
  # =========================
  det_rate <- rowMeans(detected_fun(mat))
  
  df_det <- data.frame(
    Protein = prot,
    det_rate = det_rate,
    stringsAsFactors = FALSE
  )
  
  # =========================
  # 3) 排序并生成 rank
  # =========================
  df_det <- df_det |>
    dplyr::arrange(dplyr::desc(det_rate)) |>
    dplyr::mutate(rank = dplyr::row_number())
  
  # =========================
  # 4) 分箱（仅当 color_by_bin=TRUE 时需要）
  # =========================
  if (color_by_bin) {
    if (length(bin_breaks) != 5) stop("bin_breaks 需要长度=5，例如 c(1,0.9,0.8,0.7,0.5)")
    
    df_det <- df_det |>
      dplyr::mutate(bin = dplyr::case_when(
        det_rate >= bin_breaks[1] ~ bin_lv[1],
        det_rate >= bin_breaks[2] ~ bin_lv[2],
        det_rate >= bin_breaks[3] ~ bin_lv[3],
        det_rate >= bin_breaks[4] ~ bin_lv[4],
        det_rate >= bin_breaks[5] ~ bin_lv[5],
        TRUE                      ~ bin_lv[6]
      )) |>
      dplyr::mutate(bin = factor(bin, levels = bin_lv))
    
    if (is.null(bin_colors)) {
      bin_colors <- c(
        ">=100%" = "#1b9e77",
        ">=90%"  = "#66a61e",
        ">=80%"  = "red",
        ">=70%"  = "#e6ab02",
        ">=50%"  = "skyblue",
        "<50%"   = "#7570b3"
      )
    }
    if (is.null(names(bin_colors))) {
      stop("bin_colors 建议传入带名字向量，例如 c('>=100%'='#1b9e77', ...)")
    }
  }
  
  # =========================
  # 5) x轴 breaks：每隔 x_break_by 显示一次
  # =========================
  max_rank <- max(df_det$rank, na.rm = TRUE)
  x_breaks <- seq(1, max_rank, by = x_break_by)
  
  # =========================
  # 6) 注释文本（Detected proteins ≥xx%）
  # =========================
  ann_text <- NULL
  if (show_annotation) {
    cnt <- sapply(ann_thr, function(t) sum(det_rate >= t, na.rm = TRUE))
    ann_text <- paste0(
      "Detected proteins:\n",
      paste(sprintf("≥%d%% : %d", ann_thr * 100, cnt), collapse = "\n")
    )
  }
  
  # =========================
  # 7) ✅ 标记数据准备（新增）
  # =========================
  df_mark <- NULL
  if (!is.null(mark_proteins) && length(mark_proteins) > 0) {
    df_mark <- df_det |>
      dplyr::filter(.data$Protein %in% mark_proteins) |>
      dplyr::mutate(
        y_end = .data$det_rate + mark_line_len,
        y_lab = .data$y_end
      )
    
    # 统一偏移 + 单基因偏移
    df_mark <- df_mark |>
      dplyr::mutate(x_lab = .data$rank + mark_x_nudge)
    
    if (!is.null(mark_x_nudge_map)) {
      if (is.null(names(mark_x_nudge_map))) {
        stop("mark_x_nudge_map 必须是命名向量，例如 c('ALB'=80, 'TP53'=-20)")
      }
      add_nudge <- unname(mark_x_nudge_map[df_mark$Protein])
      add_nudge[is.na(add_nudge)] <- 0
      df_mark$x_lab <- df_mark$x_lab + add_nudge
    }
    
    # 若需要，自动把y上限放宽，以容纳标签/虚线
    if (auto_expand_y_for_labels) {
      y_top_need <- max(df_mark$y_end, na.rm = TRUE)
      if (is.finite(y_top_need) && y_top_need > y_limits[2]) {
        y_limits[2] <- y_top_need + y_pad_frac
      }
    }
  }
  
  # =========================
  # 8) 绘图（柱子）
  # =========================
  p <- ggplot2::ggplot(df_det, ggplot2::aes(x = rank, y = det_rate)) +
    ggplot2::geom_col(
      ggplot2::aes(fill = if (color_by_bin) bin else NULL),
      color = NA,
      linewidth = 0.2,
      alpha = alpha,
      width = bar_width
    ) +
    ggplot2::scale_y_continuous(
      limits = y_limits,
      breaks = y_breaks,
      labels = scales::percent_format(accuracy = y_label_accuracy),
      minor_breaks = NULL
    ) +
    ggplot2::scale_x_continuous(
      breaks = x_breaks,
      labels = x_breaks,
      expand = c(0, 0)
    ) +
    ggplot2::labs(x = x_lab, y = y_lab, fill = NULL) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = if (transparent_bg) "transparent" else "white", colour = NA),
      plot.background  = ggplot2::element_rect(fill = if (transparent_bg) "transparent" else "white", colour = NA),
      legend.background = ggplot2::element_rect(fill = if (transparent_bg) "transparent" else "white", colour = NA),
      legend.box.background = ggplot2::element_rect(fill = if (transparent_bg) "transparent" else "white", colour = NA),
      
      # 坐标轴线/刻度线都删掉，但保留文字
      axis.line  = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      
      axis.text.x = ggplot2::element_text(
        color = "black",
        angle = x_angle,
        vjust = x_vjust,
        hjust = x_hjust,
        margin = ggplot2::margin(t = x_text_margin_t)
      ),
      axis.text.y = ggplot2::element_text(color = "black"),
      
      plot.margin = plot_margin,
      legend.position = if (show_legend && color_by_bin) "top" else "none"
    ) +
    ggplot2::coord_cartesian(clip = "off")
  
  # 分箱配色 / 单色
  if (color_by_bin) {
    p <- p + ggplot2::scale_fill_manual(breaks = bin_lv, values = bin_colors)
  } else {
    p <- p + ggplot2::scale_fill_identity()
    # 直接覆盖柱子统一颜色（更稳定）
    p <- ggplot2::ggplot(df_det, ggplot2::aes(x = rank, y = det_rate)) +
      ggplot2::geom_col(fill = single_color, alpha = alpha, width = bar_width) +
      ggplot2::scale_y_continuous(
        limits = y_limits, breaks = y_breaks,
        labels = scales::percent_format(accuracy = y_label_accuracy),
        minor_breaks = NULL
      ) +
      ggplot2::scale_x_continuous(breaks = x_breaks, labels = x_breaks, expand = c(0, 0)) +
      ggplot2::labs(x = x_lab, y = y_lab) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        panel.border = ggplot2::element_blank(),
        panel.background = ggplot2::element_rect(fill = if (transparent_bg) "transparent" else "white", colour = NA),
        plot.background  = ggplot2::element_rect(fill = if (transparent_bg) "transparent" else "white", colour = NA),
        
        axis.line  = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.title = ggplot2::element_blank(),
        
        axis.text.x = ggplot2::element_text(
          color = "black",
          angle = x_angle,
          vjust = x_vjust,
          hjust = x_hjust,
          margin = ggplot2::margin(t = x_text_margin_t)
        ),
        axis.text.y = ggplot2::element_text(color = "black"),
        
        plot.margin = plot_margin,
        legend.position = "none"
      ) +
      ggplot2::coord_cartesian(clip = "off")
  }
  
  # =========================
  # 9) ✅ 添加标记层（新增）
  # =========================
  if (!is.null(df_mark) && nrow(df_mark) > 0) {
    # 红点（柱顶）
    p <- p + ggplot2::geom_point(
      data = df_mark,
      ggplot2::aes(x = rank, y = det_rate),
      inherit.aes = FALSE,
      size = mark_point_size,
      color = mark_color
    ) +
      # 虚线（柱顶 -> y_end）
      ggplot2::geom_segment(
        data = df_mark,
        ggplot2::aes(x = rank, xend = rank, y = det_rate, yend = y_end),
        inherit.aes = FALSE,
        linetype = mark_line_type,
        linewidth = mark_line_width,
        color = mark_line_color
      ) +
      # 标签（虚线终点）
      ggplot2::geom_text(
        data = df_mark,
        ggplot2::aes(
          x = x_lab, y = y_lab,
          label = paste0(
            scales::percent(det_rate, accuracy = 1),  # 例如 15%
            ": ",
            Protein
          )
        ),
        inherit.aes = FALSE,
        vjust = -mark_text_vjust,
        
        size = mark_label_size,
        color = mark_label_color
      )
    
  }
  
  # =========================
  # 10) 注释框
  # =========================
  if (show_annotation && !is.null(ann_text)) {
    p <- p + ggplot2::annotate(
      "label",
      x = Inf, y = Inf,
      label = ann_text,
      hjust = ann_hjust, vjust = ann_vjust,
      size = ann_size,
      label.size = ann_label_size
    )
  }
  
  return(list(plot = p, data = df_det, mark_data = df_mark))
}


if(F){
  #1）按 bin 上色 + 标记某些蛋白
  labels <- c("ALB","SEM6B","LIRA4","PLA2R","PPIE","NENF","XPP1","HPLN3","TXTP","TM1L2","EF1D")
  
  res <- plot_detection_rate_rank(
    matrix_all,
    color_by_bin = TRUE,
    mark_proteins = labels,
    mark_line_len = 0.05,          # 5% 的虚线长度（det_rate 0~1）
    mark_x_nudge_map = c(ALB = 80) # 只对ALB标签往右挪80个rank（你自己调）
  )
  
  res$plot
  
  #2）不分箱：统一颜色 + 也照样可以标记
  res <- plot_detection_rate_rank(
    matrix_all,
    color_by_bin = FALSE,
    single_color = "#1496D4",
    mark_proteins = labels,
    mark_line_len = 0.05
  )
  
  res$plot
  
}