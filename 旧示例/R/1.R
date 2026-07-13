#### 准备 ####
setwd("C:/Work/SH/PACS/望京医院数据/7.8找HELP重叠情况/R/")
library(dplyr)
library(stringr)
library(tibble)
library(readxl)
library(tidyr)
library(ggplot2)
library(showtext)
library(patchwork)
library(factoextra)
library(FactoMineR)
library(ggrepel)
source("C:/Work/SH/code/source/CV作图.R")
source("C:/Work/SH/code/source/boxplot.R")

dir.create("../Results/", recursive = 1, showWarnings = F)
cols_group <- c(
  "1" = "#5A8BF9",
  "2" = "#39D88C",
  "3" = "#FFD400",
  "4" = "#606060"
)
cols_group3 <- c(
  "1" = "#5A8BF9",
  "2" = "#39D88C",
  "3" = "#FFD400",
  "4" = "#606060",
  "5" = "#E5843C"
)

anno_cols <- c(
  "Protein.Group",
  "Protein.Ids",
  "Protein.Names",
  "Genes",
  "First.Protein.Description"
)
############# 导入处理数据 ################
load("./整理好的数据-可以直接拿来用.RData")
original_sample_heavy <- sample_heavy
original_qc_heavy <- qc_heavy_all
original_internal_heavy <- internal_heavy_all

co_help <- rownames(original_sample_heavy)
no_overlap_help <- c(
  "LAVYQAGAR",
  "DGYLFQLLR",
  "NALALFVLPK",
  "NYNLVESLK",
  "VNHVTLSQPK"
)

intersect(co_help, no_overlap_help)

sample_heavy <- original_sample_heavy %>%
  filter(!rownames(.) %in% no_overlap_help)
qc_heavy_all <- original_qc_heavy %>% filter(!rownames(.) %in% no_overlap_help)
internal_heavy_all <- original_internal_heavy %>%
  filter(!rownames(.) %in% no_overlap_help)


low_iden_sample <- readRDS("./低鉴定量样本.RDS")
outlier_sample <- c("WJ-1-43", "WJ-1-67", "WJ-1-59", "WJ-1-83")
intersect(low_iden_sample, metadata$sample_id) # 26个
intersect(unique(c(low_iden_sample, outlier_sample)), metadata$sample_id) # 29个
## 一共剔除 29 + 9 = 38 个样本

remove_sample <- intersect(
  unique(c(low_iden_sample, outlier_sample)),
  metadata$sample_id
)

#saveRDS(remove_sample,"./第四版剔除的样本（再加9例）.RDS")

sample_heavy <- sample_heavy %>% select(-all_of(remove_sample))
metadata3 <- metadata3 %>% filter(!.$sample_id %in% remove_sample)
matrix2 <- matrix2 %>% select(-all_of(remove_sample))

if (T) {
  rm(align_samples)
  source("C:/Work/SH/code/source/对齐样本.R")
  out <- align_samples(
    exp_data = matrix2,
    group_data = metadata3,
    sample_col = "sample_id",
    group_col = c("group1"),
    #keep_groups  = list(Study = c("S1","S2")),
    group_levels = list(
      #Study = c("S1","S2"),
      #group = c("HC","ST")   # 你想要的组顺序
    ),
    sort_within = "metadata",
    verbose = TRUE
  )
}


sample_heavy_step0_boxplot <- plot_expr_boxplot(
  sample_heavy,
  metadata = metadata3,
  sample_id_col = "sample_id",
  color_by = "group1",
  order_by = "color_by",
  show_x_text = FALSE,
  #cols = cols_group4,
  show_points = F,
  log2_transform = F,
  # ymax = 20,
  title = "Raw study sample heavy \nlabeled peptides intensity"
)
ggsave(
  "../Results/sample/step0/sample-PICS分布.png",
  sample_heavy_step0_boxplot,
  width = 20,
  height = 6
)
1


source("C:/Work/SH/code/source/质谱样本名处理.R")

intersect(
  qc_internal_metadata_all$sample_id,
  names(qc_expr_all)[6:ncol(qc_expr_all)]
)

qc_expr_final <- qc_expr_final %>%
  dplyr::select(all_of(qc_metadata_all$sample_id))
qc_expr_step0 <- qc_expr_final %>% na.omit() %>% log2()


########################## 检查HELP #############
source("C:/Work/SH/somalogic/lianchuan/肽段cv 箱线图.R")
help_sample_box <- plot_peptide_cv_box(
  sample_heavy,
  value_col = "intensity",
  ylab = "log2 Heavy intensity",
  out_file = "../Results/sample/sample HELP boxplot.png"
)
help_sample_box$plot

help_qc_box <- plot_peptide_cv_box(
  qc_heavy_all,
  value_col = "intensity",
  ylab = "log2 Heavy intensity",
  out_file = "../Results/sample/qc HELP boxplot.png"
)
help_qc_box$plot

help_internal_box <- plot_peptide_cv_box(
  internal_heavy_all,
  value_col = "intensity",
  ylab = "log2 Heavy intensity",
  out_file = "../Results/sample/internal HELP boxplot.png"
)
help_internal_box$plot

############### 蛋白检出率 ##################
source("C:/Work/SH/code/source/蛋白检出率.R")
ident_rank <- plot_detection_rate_rank(
  matrix2,
  x_break_by = 200,
  y_breaks = seq(0, 1, 0.2)
)

ident_rank$plot
ggsave("../Results/蛋白检出率.png", width = 40, height = 10)
######################## step 1 #######################
############### Internal ##############
internal_heavy_step0 <- internal_heavy_all %>%
  log2() %>%
  filter(rownames(.) %in% rownames(sample_heavy)) %>%
  na.omit()
run_all <- internal_metadata_all

rm(plot_expr_boxplot)
source("C:/Work/SH/code/source/boxplot.R")
internal_h_step0_boxplot <- plot_expr_boxplot(
  internal_heavy_step0,
  metadata = run_all,
  sample_id_col = "sample_id",
  color_by = "Run",
  order_by = "color_by",
  show_x_text = FALSE,
  cols = cols_group,
  show_points = F,
  log2_transform = F,
  xlab = "Internal Calibrators HELP",
  # ymax = 20,
  title = "Raw Internal Calibrators HELP"
)
internal_h_step0_boxplot
ggsave(
  "../Results/sample/Internal-QC heavy step0-box.png",
  internal_h_step0_boxplot,
  width = 6,
  height = 5
)


#### Internal-QC HEAVY 做median normalization ####
#------------- median normalization ----------------
internal_mat1 <- as.matrix(internal_heavy_step0)

# 1) 计算每列的中位数（忽略 NA）
col_meds_internal_HELP <- apply(internal_mat1, 2, median, na.rm = TRUE)

# 2) 选定目标中位数：所有列中位数的中位数（稳健）
target_median_internal_HELP <- median(col_meds_internal_HELP, na.rm = TRUE)
target_median_internal_HELP
# 3) 每列需要加的 offset
offset_mednorm_internal_HELP <- target_median_internal_HELP -
  col_meds_internal_HELP
offset_mednorm_internal_HELP

#-------------- 计算 Intensity median --------------
# 计算每个样本的 Intensity median ##
internal_heavy_median <- data.frame(
  Sample_id = names(col_meds_internal_HELP),
  Median_intensity = col_meds_internal_HELP,
  row.names = NULL
)

# 计算每个样本的系数
internal_heavy_median$k_internal <- offset_mednorm_internal_HELP

#-------------- 矫正Internal-QC HEAVY  --------------
# “样本名 -> k_sam”的命名向量
k_vec_internal_heavy <- setNames(
  internal_heavy_median$k_internal,
  internal_heavy_median$Sample_id
)

# 按列乘（sweep 会自动对齐列）
internal_heavy_step1 <- sweep(
  internal_heavy_step0,
  2,
  k_vec_internal_heavy,
  `+`
)

internal_h_step1_boxplot <- plot_expr_boxplot(
  internal_heavy_step1,
  metadata = run_all,
  sample_id_col = "sample_id",
  color_by = "Run",
  order_by = "color_by",
  show_x_text = FALSE,
  cols = cols_group,
  show_points = F,
  log2_transform = F,
  xlab = "Internal Calibrators HELP",
  # ymax = 20,
  title = "Step 1 Internal Calibrators HELP"
)
internal_h_step1_boxplot
ggsave(
  "../Results/sample/Internal-QC heavy step1-box.png",
  internal_h_step1_boxplot,
  width = 6,
  height = 5
)
############### QC ##############
#sample_heavy_bak <- sample_heavy
#sample_heavy <- sample_heavy_bak
#sele_pics2 <- readRDS("../sele_pics2.RDS")
#sample_heavy <- sample_heavy %>% filter(rownames(.) %in% sele_pics2)
qc_heavy_step0 <- qc_heavy_all %>%
  log2() %>%
  filter(rownames(.) %in% rownames(sample_heavy)) %>%
  na.omit()

rm(plot_expr_boxplot)
source("C:/Work/SH/code/source/boxplot.R")
qc_h_step0_boxplot <- plot_expr_boxplot(
  qc_heavy_step0,
  metadata = qc_metadata_all,
  sample_id_col = "sample_id",
  color_by = "Plate",
  order_by = "color_by",
  show_x_text = FALSE,
  cols = cols_group,
  show_points = F,
  log2_transform = F,
  xlab = "QC sample HELP",
  # ymax = 20,
  title = "Raw QC HELP"
)
#qc_h_step0_boxplot

ggsave(
  "../Results/sample/QC heavy step0-box.png",
  qc_h_step0_boxplot,
  width = 6,
  height = 5
)

#### QC HEAVY 做median normalization ####
#------------- median normalization ----------------
qc_mat1_HELP <- as.matrix(qc_heavy_step0)

# 1) 计算每列的中位数（忽略 NA）
col_meds_qc_HELP <- apply(qc_mat1_HELP, 2, median, na.rm = TRUE)

# 2) 选定目标中位数：所有列中位数的中位数（稳健）
target_median_qc_HELP <- median(col_meds_qc_HELP, na.rm = TRUE)
target_median_qc_HELP
# 3) 每列需要加的 offset
offset_mednorm_qc_HELP <- target_median_qc_HELP - col_meds_qc_HELP
offset_mednorm_qc_HELP

#-------------- 计算  QC HEAVY median --------------
# 计算每个样本的 Intensity median ##
qc_heavy_median <- data.frame(
  Sample_id = names(col_meds_qc_HELP),
  Median_intensity = col_meds_qc_HELP,
  row.names = NULL
)

# 计算每个样本的系数
qc_heavy_median$k_qc <- offset_mednorm_qc_HELP


#-------------- 矫正 QC HEAVY --------------
# “样本名 -> k_sam”的命名向量
k_vec_qc_heavy <- setNames(qc_heavy_median$k_qc, qc_heavy_median$Sample_id)

# 按列乘（sweep 会自动对齐列）
qc_heavy_step1 <- sweep(qc_heavy_step0, 2, k_vec_qc_heavy, `+`)


qc_h_step1_boxplot <- plot_expr_boxplot(
  qc_heavy_step1,
  metadata = qc_metadata_all,
  sample_id_col = "sample_id",
  color_by = "Plate",
  order_by = "color_by",
  show_x_text = FALSE,
  cols = cols_group3,
  show_points = F,
  log2_transform = F,
  xlab = "QC sample HELP",
  # ymax = 20,
  title = "Step 2.1 QC HELP"
)
qc_h_step1_boxplot
ggsave(
  "../Results/sample/QC heavy step1-box.png",
  qc_h_step1_boxplot,
  width = 6,
  height = 5
)

#-------------- 先对 QC 进行 step 1 ------------
internal_heavy_median
internal_heavy_median$run <- run_all$Run

median_k_internal_by_run <- internal_heavy_median %>%
  group_by(run) %>%
  summarise(
    n_qc = n(), # 每个 run 有多少个 QC
    median_k_internal = median(k_internal, na.rm = TRUE), # 你要的：组内 k_internal 中位数
    median_intensity_run = median(Median_intensity, na.rm = TRUE), # 可选：组内强度中位数
    .groups = "drop"
  ) %>%
  arrange(run)

median_k_internal_by_run


run2k <- median_k_internal_by_run %>%
  transmute(run = as.character(run), k = as.numeric(median_k_internal)) %>%
  tibble::deframe() # 得到命名向量：names=run, values=k

qc_expr_step0 <- qc_expr_final %>% na.omit() %>% log2()

meta_use_qc <- qc_metadata_all %>%
  mutate(sample_id = as.character(sample_id), Run = as.character(Run)) %>%
  filter(sample_id %in% colnames(qc_expr_step0)) %>%
  distinct(sample_id, Run)

# 每个 sample 对应的 k（命名向量：names=sample_id）
qc2k_run <- run2k[meta_use_qc$Run]
names(qc2k_run) <- meta_use_qc$sample_id

# 按矩阵列顺序对齐
k_vec_qc_run <- qc2k_run[colnames(qc_expr_step0)]

if (any(is.na(k_vec_qc_run))) {
  warning(
    "这些样本没有匹配到 Run 的 median_k_internal：\n",
    paste(names(k_vec_qc_run)[is.na(k_vec_qc_run)], collapse = ", ")
  )
}
k_vec_qc_run
qc_expr_step1 <- sweep(qc_expr_step0, 2, k_vec_qc_run, FUN = "+")
#-------------- 用 PICS 矫正 QC protein  --------------
qc_heavy_median$plate <- qc_metadata_all$Plate
median_k_qc_by_plate_for_qc <- qc_heavy_median %>%
  group_by(plate) %>%
  summarise(
    n_qc = n(), # 每个 plate 有多少个 QC
    median_k_qc = median(k_qc, na.rm = TRUE), # 你要的：组内 k_qc 中位数
    median_intensity_plate = median(Median_intensity, na.rm = TRUE), # 可选：组内强度中位数
    .groups = "drop"
  ) %>%
  arrange(plate)

plate2k_qc <- median_k_qc_by_plate_for_qc %>%
  transmute(plate = as.character(plate), k = as.numeric(median_k_qc)) %>%
  tibble::deframe() # 得到命名向量：names=run, values=k
meta_use_qc <- qc_metadata_all %>%
  mutate(sample_id = as.character(sample_id), Plate = as.character(Plate)) %>%
  filter(sample_id %in% colnames(qc_expr_final)) %>%
  distinct(sample_id, Plate)

# 每个 sample 对应的 k（命名向量：names=sample_id）
qc2k <- plate2k_qc[meta_use_qc$Plate]
names(qc2k) <- meta_use_qc$sample_id

# 按矩阵列顺序对齐
k_vec_qc <- qc2k[colnames(qc_expr_final)]

if (any(is.na(k_vec_qc))) {
  warning(
    "这些样本没有匹配到 Plate 的 median_k_qc_by_plate_for_qc：\n",
    paste(names(k_vec_qc)[is.na(k_vec_qc)], collapse = ", ")
  )
}
k_vec_qc
qc_expr_step2 <- sweep(qc_expr_step1, 2, k_vec_qc, FUN = "+")


#------------- 对 QC protein 进行 median normalization ----------------
qc_mat1 <- as.matrix(qc_expr_step2)

# 1) 计算每列的中位数（忽略 NA）
col_meds_qc <- apply(qc_mat1, 2, median, na.rm = TRUE)

# 2) 选定目标中位数：所有列中位数的中位数（稳健）
target_median_qc <- median(col_meds_qc, na.rm = TRUE)
target_median_qc

# 3) 每列需要加的 offset
offset_mednorm_qc <- target_median_qc - col_meds_qc
offset_mednorm_qc

#-------------- 计算 Intensity median --------------
# 计算每个样本的 Intensity median ##
qc_intensity_median <- data.frame(
  Sample_id = names(col_meds_qc),
  Median_intensity = col_meds_qc,
  row.names = NULL
)
qc_intensity_median
# 计算每个样本的系数
qc_intensity_median$k_qc <- offset_mednorm_qc
#qc_intensity_median_no21 <- qc_intensity_median
#-------------- 矫正 QC protein  --------------
# “样本名 -> k_sam”的命名向量
k_vec_qc <- setNames(qc_intensity_median$k_qc, qc_intensity_median$Sample_id)

# 按列乘（sweep 会自动对齐列）
qc_expr_step22 <- sweep(qc_expr_step2, 2, k_vec_qc, `+`)
#qc_expr_step22_no21 <- sweep(qc_expr_step1, 2, k_vec_qc, `+`)

qc_expr_step1_boxplot <- plot_expr_boxplot(
  qc_expr_step1,
  metadata = qc_metadata_all,
  sample_id_col = "sample_id",
  color_by = "Plate",
  order_by = "color_by",
  show_x_text = FALSE,
  cols = cols_group3,
  show_points = F,
  log2_transform = F,
  xlab = "QC sample protein",
  # ymax = 20,
  title = "Step 1 QC protein intensity"
)
qc_expr_step1_boxplot
ggsave(
  "../Results/sample/QC protein step1-box.png",
  qc_expr_step1_boxplot,
  width = 6,
  height = 5
)

qc_expr_step0_boxplot <- plot_expr_boxplot(
  qc_expr_final %>% na.omit() %>% log2(),
  metadata = qc_metadata_all,
  sample_id_col = "sample_id",
  color_by = "Plate",
  order_by = "color_by",
  show_x_text = FALSE,
  cols = cols_group3,
  show_points = F,
  log2_transform = F,
  xlab = "QC sample protein",
  # ymax = 20,
  title = "Raw QC protein intensity"
)
qc_expr_step0_boxplot
ggsave(
  "../Results/sample/QC protein step0-box.png",
  qc_expr_step0_boxplot,
  width = 6,
  height = 5
)

qc_expr_step2_boxplot <- plot_expr_boxplot(
  qc_expr_step2,
  metadata = qc_metadata_all,
  sample_id_col = "sample_id",
  color_by = "Plate",
  order_by = "color_by",
  show_x_text = FALSE,
  cols = cols_group3,
  show_points = F,
  log2_transform = F,
  xlab = "QC sample protein",
  # ymax = 20,
  title = "Step 2.1 QC protein intensity"
)
qc_expr_step2_boxplot
ggsave(
  "../Results/sample/QC protein step2-box.png",
  qc_expr_step2_boxplot,
  width = 6,
  height = 5
)

qc_expr_step22_boxplot <- plot_expr_boxplot(
  qc_expr_step22,
  metadata = qc_metadata_all,
  sample_id_col = "sample_id",
  color_by = "Plate",
  order_by = "color_by",
  show_x_text = FALSE,
  cols = cols_group3,
  show_points = F,
  log2_transform = F,
  xlab = "QC sample protein",
  # ymax = 20,
  title = "Step 2.2 QC protein intensity"
)
qc_expr_step22_boxplot
ggsave(
  "../Results/sample/QC protein step22-box.png",
  qc_expr_step22_boxplot,
  width = 6,
  height = 5
)
#-------------- QC heavy PICS & protein 的CV --------------------
rm(plot_cv_multi_stage)
source("C:/Work/SH/somalogic/lianchuan/多组cv比较_script.R")
cv_qc_all <- plot_cv_multi_stage(
  mat_list = list(
    `QC raw HELP` = qc_heavy_step0,
    `QC HELP step 2.1` = qc_heavy_step1,
    `QC raw protein` = qc_expr_step0,
    `QC protein step 1` = qc_expr_step1,
    `QC protein step 2.1` = qc_expr_step2,
    `QC protein step 2.2` = qc_expr_step22
  ),
  metadata = qc_metadata_all,
  group_col = "Plate",
  cols_stage = c(
    `QC raw HELP` = "#5A8BF9",
    `QC HELP step 2.1` = "#39D88C",
    `QC raw protein` = "#FFD400",
    `QC protein step 1` = "#606060",
    `QC protein step 2.1` = "#E5843C",
    `QC protein step 2.2` = "#CCE007"
  ),
  cv_cap = NULL,
  #title = "Step 0 vs Step 1 vs Step 2 vs Step 3.1 vs Step 3.2",
  out_file = "../Results/sample/QC-CV(HELP & protein).png",
  width = 20,
  height = 5
)
cv_qc_all$cv_all
#cv_qc_all$plot

############### Sample ##############
sample_heavy_step0 <- sample_heavy

rm(plot_expr_boxplot)
source("C:/Work/SH/code/source/boxplot.R")
sample_h_step0_boxplot <- plot_expr_boxplot(
  sample_heavy_step0,
  metadata = metadata3,
  sample_id_col = "sample_id",
  color_by = "group1",
  order_by = "color_by",
  show_x_text = T,
  cols = cols_group,
  show_points = F,
  log2_transform = F,
  xlab = "Study Sample HELP",
  # ymax = 20,
  title = "Raw sample HELP"
)
#sample_h_step0_boxplot
ggsave(
  "../Results/sample/sample heavy step0-box—tmp.png",
  sample_h_step0_boxplot,
  width = 12,
  height = 5
)


get_group_offsets <- function(
  col_meds,
  metadata,
  sample_col = 1,
  group_col = "group",
  mat_cols = NULL
) {
  # 提取样本名和分组
  sample_ids <- as.character(metadata[[sample_col]])
  group_ids <- as.character(metadata[[group_col]])

  # 每个样本对应的列中位数
  sample_meds <- col_meds[sample_ids]

  # 计算每组 target median
  group_target <- tapply(sample_meds, group_ids, median, na.rm = TRUE)

  # 每个样本的 offset
  offsets <- group_target[group_ids] - sample_meds
  names(offsets) <- sample_ids

  # 如果指定 mat_cols，则按指定顺序输出
  if (!is.null(mat_cols)) {
    offsets <- offsets[mat_cols]
  }

  return(offsets)
}
#### sample HEAVY 做median normalization ####
#-------------- 计算 Intensity median --------------
mat1 <- as.matrix(sample_heavy_step0)

# 1) 计算每列的中位数（忽略 NA）
col_meds_sam_HELP <- apply(mat1, 2, median, na.rm = TRUE)
offset_mednorm_sam_HELP <- get_group_offsets(
  col_meds = col_meds_sam_HELP,
  metadata = metadata3,
  sample_col = 1,
  group_col = "group1",
  mat_cols = colnames(mat1)
)

offset_mednorm_sam_HELP
# 计算每个样本的 Intensity median ##
sample_heavy_median <- data.frame(
  Sample_id = names(col_meds_sam_HELP),
  Median_intensity = col_meds_sam_HELP,
  row.names = NULL
)

# 计算每个样本的系数
sample_heavy_median$k_sample <- offset_mednorm_sam_HELP

#-------------- 矫正SAMPLE HEAVY  --------------
# “样本名 -> k_sam”的命名向量
k_vec_sample_heavy <- setNames(
  sample_heavy_median$k_sample,
  sample_heavy_median$Sample_id
)

# 按列乘（sweep 会自动对齐列）
sample_heavy_step1 <- sweep(sample_heavy_step0, 2, k_vec_sample_heavy, `+`)

sample_h_step1_boxplot <- plot_expr_boxplot(
  sample_heavy_step1,
  metadata = metadata3,
  sample_id_col = "sample_id",
  color_by = "group1",
  order_by = "color_by",
  show_x_text = FALSE,
  cols = cols_group,
  show_points = F,
  log2_transform = F,
  xlab = "Study Sample HELP",
  # ymax = 20,
  title = "Step 3.1 sample HELP"
)
sample_h_step1_boxplot
ggsave(
  "../Results/sample/sample heavy step1-box.png",
  sample_h_step1_boxplot,
  width = 12,
  height = 5
)

#################### step 2 ################
#### 用上面的数据给sample做矫正 ###
#-------------- 计算 --------------
qc_intensity_median
qc_intensity_median$plate <- qc_metadata_all$Plate
#qc_intensity_median_no21$plate <- qc_metadata_all$Plate

median_k_internal_by_run

median_k_qc_by_plate <- qc_intensity_median %>%
  group_by(plate) %>%
  summarise(
    n_qc = n(), # 每个 plate 有多少个 QC
    median_k_qc = median(k_qc, na.rm = TRUE), # 你要的：组内 k_qc 中位数
    median_intensity_plate = median(Median_intensity, na.rm = TRUE), # 可选：组内强度中位数
    .groups = "drop"
  ) %>%
  arrange(plate)

median_k_qc_by_plate

#median_k_qc_by_plate2 <- qc_intensity_median_no21 %>%
#  group_by(plate) %>%
#  summarise(
#    n_qc = n(),  # 每个 plate 有多少个 QC
#    median_k_qc = median(k_qc, na.rm = TRUE),  # 你要的：组内 k_qc 中位数
#    median_intensity_plate = median(Median_intensity, na.rm = TRUE),  # 可选：组内强度中位数
#    .groups = "drop"
#  ) %>%
#  arrange(plate)
#
#median_k_qc_by_plate2

#---------- 开始矫正 -----------
median_k_internal_by_run
median_k_qc_by_plate
#median_k_qc_by_plate2

sample_intensity_step0 <- log2(matrix2)
metadata3

#------------- RUN ----------------

meta_use <- metadata3 %>%
  mutate(sample_id = as.character(sample_id), Run = as.character(Run)) %>%
  filter(sample_id %in% colnames(sample_intensity_step0)) %>%
  distinct(sample_id, Run)

# 每个 sample 对应的 k（命名向量：names=sample_id）
sample2k <- run2k[meta_use$Run]
names(sample2k) <- meta_use$sample_id

# 按矩阵列顺序对齐
k_vec <- sample2k[colnames(sample_intensity_step0)]

if (any(is.na(k_vec))) {
  warning(
    "这些样本没有匹配到 Run 的 median_k_internal：\n",
    paste(names(k_vec)[is.na(k_vec)], collapse = ", ")
  )
}

sample_intensity_step1 <- sweep(sample_intensity_step0, 2, k_vec, FUN = "+")

# sample_intensity_step1 就是校正后的矩阵

#------------- Plate ----------------
plate2k <- median_k_qc_by_plate %>% # median_k_qc_by_plate2
  transmute(plate = as.character(plate), k = as.numeric(median_k_qc)) %>%
  tibble::deframe() # 命名向量：names=plate, values=k

meta_use2 <- metadata3 %>%
  mutate(sample_id = as.character(sample_id), Plate = as.character(Plate)) %>%
  filter(sample_id %in% colnames(sample_intensity_step1)) %>%
  distinct(sample_id, Plate)

# 每个 sample 对应的 plate k（命名向量：names=sample_id）
sample2k_plate <- plate2k[meta_use2$Plate]
names(sample2k_plate) <- meta_use2$sample_id

# 按矩阵列顺序对齐
k_plate_vec <- sample2k_plate[colnames(sample_intensity_step1)]

if (any(is.na(k_plate_vec))) {
  warning(
    "这些样本没有匹配到 Plate 的 median_k_qc：\n",
    paste(names(k_plate_vec)[is.na(k_plate_vec)], collapse = ", ")
  )
}

sample_intensity_step2 <- sweep(
  sample_intensity_step1,
  2,
  k_plate_vec,
  FUN = "+"
)
#sample_intensity_step2_no21 <- sweep(sample_intensity_step1, 2, k_plate_vec, FUN = "+")

# sample_intensity_step2 即：先按 Run 校正，再按 Plate 校正后的矩阵

#------------- Sample ----------------
sample2k_sample <- sample_heavy_median %>%
  transmute(
    sample_id = as.character(Sample_id),
    k_sample = as.numeric(k_sample)
  ) %>%
  tibble::deframe() # 命名向量：names=sample_id, values=k_sample

k_sample_vec <- sample2k_sample[colnames(sample_intensity_step2)]

if (any(is.na(k_sample_vec))) {
  warning(
    "这些样本在 sample_heavy_median 里找不到 k_sample：\n",
    paste(
      colnames(sample_intensity_step2)[is.na(k_sample_vec)],
      collapse = ", "
    )
  )
}

sample_intensity_step3 <- sweep(
  sample_intensity_step2,
  2,
  k_sample_vec,
  FUN = "+"
)
#sample_intensity_step3_no21 <- sweep(sample_intensity_step2_no21, 2, k_sample_vec, FUN = "+")

# sample_intensity_step3 就是最终矩阵

#------------- median normalization ----------------
mat3 <- as.matrix(sample_intensity_step3)
#mat3 <- as.matrix(sample_intensity_step3_no21)

# 1) 计算每列的中位数（忽略 NA）
col_meds <- apply(mat3, 2, median, na.rm = TRUE)

offset_mednorm <- get_group_offsets(
  col_meds = col_meds,
  metadata = metadata3,
  sample_col = 1,
  group_col = "group1",
  mat_cols = colnames(mat3)
)

offset_mednorm


# 计算每个样本的 Intensity median ##
sample_median <- data.frame(
  Sample_id = names(col_meds),
  Median_intensity = col_meds,
  row.names = NULL
)

# 计算每个样本的系数
sample_median$k_sample <- offset_mednorm

# 4) 对每列整体加 offset
sample_intensity_step4 <- sweep(mat3, 2, offset_mednorm, FUN = "+") %>%
  as.data.frame()
#sample_intensity_step4_no21 <- sweep(mat3, 2, offset_mednorm, FUN = "+") %>%
#  as.data.frame()

#---------------- 可视化校准后的 sample Intensity -----------------
### CV ###
### CV ###
### CV ###
rm(plot_cv_by_metadata)
source("C:/Work/SH/code/source/CV作图.R")
cv_0 <- plot_cv_by_metadata(
  mat = sample_intensity_step0,
  metadata = metadata3,
  sample_id_col = "sample_id",
  group_col = "group1",
  #cols = cols_group2,
  min_detect_rate = 0.70, # ★组内>=70%非NA才算
  min_non_na = 3, # 同时至少>=3个非NA
  #cv_cap = 200,             # 可选：避免少数极端CV拉高y轴
  title = "Raw sample intensity",
  out_prefix = "../Results/sample/step0/CV", # ★只给前缀
  save_formats = c("png"),
  width = 8,
  height = 6,
  dpi = 300
)
cv_1 <- plot_cv_by_metadata(
  mat = sample_intensity_step1,
  metadata = metadata3,
  sample_id_col = "sample_id",
  group_col = "group1",
  #cols = cols_group2,
  min_detect_rate = 0.70, # ★组内>=70%非NA才算
  min_non_na = 3, # 同时至少>=3个非NA
  #cv_cap = 200,             # 可选：避免少数极端CV拉高y轴
  title = "Step 1 sample intensity",
  out_prefix = "../Results/sample/step1/CV", # ★只给前缀
  save_formats = c("png"),
  width = 8,
  height = 6,
  dpi = 300
)
cv_2 <- plot_cv_by_metadata(
  mat = sample_intensity_step2,
  metadata = metadata3,
  sample_id_col = "sample_id",
  group_col = "group1",
  #cols = cols_group2,
  min_detect_rate = 0.70, # ★组内>=70%非NA才算
  min_non_na = 3, # 同时至少>=3个非NA
  #cv_cap = 200,             # 可选：避免少数极端CV拉高y轴
  title = "Step 2 sample intensity",
  out_prefix = "../Results/sample/step2/CV", # ★只给前缀
  save_formats = c("png"),
  width = 8,
  height = 6,
  dpi = 300
)
cv_3 <- plot_cv_by_metadata(
  mat = sample_intensity_step3,
  metadata = metadata3,
  sample_id_col = "sample_id",
  group_col = "group1",
  #cols = cols_group2,
  min_detect_rate = 0.70, # ★组内>=70%非NA才算
  min_non_na = 3, # 同时至少>=3个非NA
  #cv_cap = 200,             # 可选：避免少数极端CV拉高y轴
  title = "Step 3.1 sample intensity",
  out_prefix = "../Results/sample/step3/CV", # ★只给前缀
  save_formats = c("png"),
  width = 8,
  height = 6,
  dpi = 300
)
cv_4 <- plot_cv_by_metadata(
  mat = sample_intensity_step4,
  metadata = metadata3,
  sample_id_col = "sample_id",
  group_col = "group1",
  #cols = cols_group2,
  min_detect_rate = 0.70, # ★组内>=70%非NA才算
  min_non_na = 3, # 同时至少>=3个非NA
  #cv_cap = 200,             # 可选：避免少数极端CV拉高y轴
  title = "Step 3.2 sample intensity",
  out_prefix = "../Results/sample/step4/CV-single", # ★只给前缀
  save_formats = c("png"),
  width = 8,
  height = 6,
  dpi = 300
)


rm(plot_cv_multi_stage)
source("C:/Work/SH/somalogic/lianchuan/多组cv比较_script.R")
cv_sample_all <- plot_cv_multi_stage(
  mat_list = list(
    `Step 0` = sample_intensity_step0,
    `Step 1` = sample_intensity_step1,
    `Step 2` = sample_intensity_step2,
    `Step 3.1` = sample_intensity_step3,
    `Step 3.2` = sample_intensity_step4
  ),
  metadata = metadata3,
  group_col = "group1",
  cols_stage = c(
    `Step 0` = "#5A8BF9",
    `Step 1` = "#39D88C",
    `Step 2` = "#FFD400",
    `Step 3.1` = "#606060",
    `Step 3.2` = "#E5843C",
    `Step 4` = "#CCE007"
  ),
  cv_cap = NULL,
  #title = "Step 0 vs Step 1 vs Step 2 vs Step 3.1 vs Step 3.2",
  out_file = "../Results/sample/step5/CV.png",
  width = 16,
  height = 5
)
#cv_sample_all$plot

### boxplot ###
### boxplot ###
### boxplot ###
rm(plot_expr_boxplot)
source("C:/Work/SH/code/source/boxplot.R")
sam_step0_boxplot <- plot_expr_boxplot(
  sample_intensity_step0,
  metadata = metadata3,
  sample_id_col = "sample_id",
  color_by = "group1",
  order_by = "color_by",
  show_x_text = FALSE,
  #cols = cols_group1,
  show_points = F,
  log2_transform = F,
  xlab = "group1 Sample",
  fill_alpha = 0.7, # ★颜色本身透明度
  # ymax = 20,
  title = "Raw sample intensity"
)
ggsave(
  "../Results/sample/step0/boxplot.png",
  sam_step0_boxplot,
  width = 12,
  height = 5
)

sam_step1_boxplot <- plot_expr_boxplot(
  sample_intensity_step1,
  metadata = metadata3,
  sample_id_col = "sample_id",
  color_by = "group1",
  order_by = "color_by",
  show_x_text = FALSE,
  #cols = cols_group1,
  show_points = F,
  log2_transform = F,
  xlab = "group1 Sample",
  fill_alpha = 0.7, # ★颜色本身透明度
  # ymax = 20,
  title = "Step 1 sample intensity"
)

ggsave(
  "../Results/sample/step1/boxplot.png",
  sam_step1_boxplot,
  width = 12,
  height = 5
)
1
sam_step2_boxplot <- plot_expr_boxplot(
  sample_intensity_step2,
  metadata = metadata3,
  sample_id_col = "sample_id",
  color_by = "group1",
  order_by = "color_by",
  show_x_text = FALSE,
  #cols = cols_group1,
  show_points = F,
  log2_transform = F,
  xlab = "group1 Sample",
  fill_alpha = 0.7, # ★颜色本身透明度
  # ymax = 20,
  title = "Step 2 sample intensity"
)
ggsave(
  "../Results/sample/step2/boxplot.png",
  sam_step2_boxplot,
  width = 12,
  height = 5
)
1
sam_step3_boxplot <- plot_expr_boxplot(
  sample_intensity_step3,
  metadata = metadata3,
  sample_id_col = "sample_id",
  color_by = "group1",
  order_by = "color_by",
  show_x_text = FALSE,
  #cols = cols_group1,
  show_points = F,
  log2_transform = F,
  xlab = "group1 Sample",
  fill_alpha = 0.7, # ★颜色本身透明度
  # ymax = 20,
  title = "Step 3.1 sample intensity"
)
ggsave(
  "../Results/sample/step3/boxplot.png",
  sam_step3_boxplot,
  width = 12,
  height = 5
)
1
sam_step4_boxplot <- plot_expr_boxplot(
  sample_intensity_step4,
  metadata = metadata3,
  sample_id_col = "sample_id",
  color_by = "group1",
  order_by = "color_by",
  show_x_text = FALSE,
  #cols = cols_group1,
  show_points = F,
  log2_transform = F,
  xlab = "group1 Sample",
  fill_alpha = 0.7, # ★颜色本身透明度
  # ymax = 20,
  title = "Step 3.2 sample intensity"
)
ggsave(
  "../Results/sample/step4/boxplot—single.png",
  sam_step4_boxplot,
  width = 12,
  height = 5
)
1

ggsave(
  "../Results/sample/step5/箱线图.png",
  (sam_step0_boxplot +
    sam_step1_boxplot +
    sam_step2_boxplot +
    sam_step3_boxplot +
    sam_step4_boxplot),
  width = 60,
  height = 20,
  limitsize = FALSE
)

### PCA ###
### PCA ###
### PCA ###
if (F) {
  rm(plot_pca_samples)
  source("C:/Work/SH/code/source/PCA-plot.R")
  pca_step0 <- plot_pca_samples(
    matrix2, #palette = cols_group2,
    addEllipses = TRUE,
    metadata3$group1,
    show_labels = F,
    palette_alpha = 0.7,
    title = "Raw sample intensity"
  )
  #pca_step0
  ggsave("../Results/sample/step0/PCA.png", pca_step0, width = 8, height = 5)

  pca_step1 <- plot_pca_samples(
    sample_intensity_step1,
    palette = cols_group2,
    metadata3$group,
    show_labels = F,
    palette_alpha = 0.7,
    title = "Step 1 sample intensity"
  )
  #pca_step1
  ggsave("../Results/sample/step1/PCA.png", pca_step1, width = 5, height = 5)
  pca_step2 <- plot_pca_samples(
    sample_intensity_step2,
    palette = cols_group2,
    metadata3$group,
    show_labels = F,
    palette_alpha = 0.7,
    title = "Step 2 sample intensity"
  )
  ggsave("../Results/sample/step2/PCA.png", pca_step2, width = 5, height = 5)
  pca_step3 <- plot_pca_samples(
    sample_intensity_step3,
    palette = cols_group2,
    metadata3$group,
    show_labels = F,
    palette_alpha = 0.7,
    title = "Step 3 sample intensity"
  )
  ggsave("../Results/sample/step3/PCA.png", pca_step3, width = 5, height = 5)
  pca_step4 <- plot_pca_samples(
    sample_intensity_step4,
    palette = cols_group2,
    metadata3$group,
    show_labels = F,
    palette_alpha = 0.7,
    title = "Step 4 sample intensity"
  )
  ggsave(
    "../Results/sample/step4/PCA-single.png",
    pca_step4,
    width = 5,
    height = 5
  )
}

rm(plot_pca_samples)
source("C:/Work/SH/code/source/PCA-plot.R")
mat_list <- list(
  Step0 = sample_intensity_step0,
  Step1 = sample_intensity_step1,
  Step2 = sample_intensity_step2,
  Step3 = sample_intensity_step3,
  Step4 = sample_intensity_step4
)

lims <- get_common_pca_limits(mat_list, axes = c(1, 2), symmetric = TRUE)

pca_step0 <- plot_pca_samples(
  sample_intensity_step0,
  group = metadata3$group1,
  palette_alpha = 0.7,
  title = "Raw sample intensity",
  addEllipses = F,
  show_legend = FALSE,
  axis_title_style = "blank", # 或 "blank"
  xlim = lims$xlim,
  ylim = lims$ylim,
  axis_on_edge = TRUE,
  remove_origin = TRUE
)
ggsave("../Results/sample/step0/PCA.png", pca_step0, width = 5, height = 5)
pca_step1 <- plot_pca_samples(
  sample_intensity_step1,
  group = metadata3$group1,
  palette_alpha = 0.7,
  title = "Step 1 sample intensity",
  addEllipses = F,
  show_legend = FALSE,
  axis_title_style = "blank", # "PC"或 "blank"
  xlim = lims$xlim,
  ylim = lims$ylim,
  axis_on_edge = TRUE,
  remove_origin = TRUE
)
ggsave("../Results/sample/step1/PCA.png", pca_step1, width = 5, height = 5)

pca_step2 <- plot_pca_samples(
  sample_intensity_step2,
  group = metadata3$group1,
  palette_alpha = 0.7,
  title = "Step 2 sample intensity",
  addEllipses = F,
  show_legend = FALSE,
  axis_title_style = "blank", # 或 "blank"
  xlim = lims$xlim,
  ylim = lims$ylim,
  axis_on_edge = TRUE,
  remove_origin = TRUE
)
ggsave("../Results/sample/step2/PCA.png", pca_step2, width = 5, height = 5)
pca_step3 <- plot_pca_samples(
  sample_intensity_step3,
  group = metadata3$group1,
  palette_alpha = 0.7,
  title = "Step 3.1 sample intensity",
  addEllipses = F,
  show_legend = FALSE,
  axis_title_style = "blank", # 或 "blank"
  xlim = lims$xlim,
  ylim = lims$ylim,
  axis_on_edge = TRUE,
  remove_origin = TRUE
)
ggsave("../Results/sample/step3/PCA.png", pca_step3, width = 5, height = 5)

pca_step4 <- plot_pca_samples(
  sample_intensity_step4,
  group = metadata3$group1,
  palette_alpha = 0.7,
  title = "Step 3.2 sample intensity",
  addEllipses = F,
  show_legend = T,
  axis_title_style = "blank", # 或 "blank"
  xlim = lims$xlim,
  ylim = lims$ylim,
  axis_on_edge = TRUE,
  remove_origin = TRUE
)
ggsave(
  "../Results/sample/step4/PCA-single.png",
  pca_step4,
  width = 5,
  height = 5
)

pca_step0 | pca_step1 | pca_step2 | pca_step3 | pca_step4
ggsave("../Results/sample/step5/PCA.png", width = 25, height = 5)
#ggsave("../Results/sample/step4/PCA.pdf", width = 25, height = 5)

### Density ###
### Density ###
### Density ###
rm(plot_plate_density)
source("C:/Work/SH/somalogic/lianchuan/分布密度图.R")
density_step0_sample <- plot_plate_density(
  sample_intensity_step0,
  metadata3,
  plate_col = "group1",
  #cols = cols_group2,
  color_alpha = 0.3,
  mode = "values",
  linewidth = 0.3,
  title = "Raw sample density",
  n_per_plate = 200000,
  color_by = "sample_by_plate",
  legend = TRUE
)
ggsave(
  "../Results/sample/step0/Density.png",
  density_step0_sample,
  width = 10,
  height = 5
)


density_step1_sample <- plot_plate_density(
  sample_intensity_step1,
  metadata3,
  plate_col = "group1",
  #cols = cols_group2,
  color_alpha = 0.3,
  linewidth = 0.3,
  mode = "values",
  title = "Step 1 sample density",
  n_per_plate = 200000,
  color_by = "sample_by_plate",
  legend = TRUE
)
#density_step1_sample
ggsave(
  "../Results/sample/step1/Density.png",
  density_step1_sample,
  width = 10,
  height = 5
)
density_step2_sample <- plot_plate_density(
  sample_intensity_step2,
  metadata3,
  plate_col = "group1",
  #cols = cols_group2,
  color_alpha = 0.3,
  linewidth = 0.3,
  mode = "values",
  title = "Step 2 sample density",
  n_per_plate = 200000,
  color_by = "sample_by_plate",
  legend = TRUE
)
ggsave(
  "../Results/sample/step2/Density.png",
  density_step2_sample,
  width = 10,
  height = 5
)

density_step3_sample <- plot_plate_density(
  sample_intensity_step3,
  metadata3,
  plate_col = "group1",
  #cols = cols_group2,
  color_alpha = 0.3,
  linewidth = 0.3,
  mode = "values",
  title = "Step 3.1 sample density",
  n_per_plate = 200000,
  color_by = "sample_by_plate",
  legend = TRUE
)
ggsave(
  "../Results/sample/step3/Density.png",
  density_step3_sample,
  width = 10,
  height = 5
)

density_step4_sample <- plot_plate_density(
  sample_intensity_step4,
  metadata3,
  plate_col = "group1",
  #cols = cols_group2,
  color_alpha = 0.3,
  linewidth = 0.3,
  mode = "values",
  title = "Step 3.2 sample density",
  n_per_plate = 200000,
  color_by = "sample_by_plate",
  legend = TRUE
)
ggsave(
  "../Results/sample/step4/Density-single.png",
  density_step4_sample,
  width = 10,
  height = 5
)


density_step0_sample /
  density_step1_sample /
  density_step2_sample /
  density_step3_sample /
  density_step4_sample
ggsave(
  "../Results/sample/step5/密度图.png",
  (density_step0_sample /
    density_step1_sample /
    density_step2_sample /
    density_step3_sample /
    density_step4_sample),
  width = 20,
  height = 30
)
save.image("./第四版结果-剔除38个样本.RData")
