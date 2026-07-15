# PACS strategy benchmark

从仓库根目录运行：

```bash
Rscript script/run_all.R
```

依赖包：

```r
install.packages(c("data.table", "readxl", "ggplot2"))
```

当前实现的策略：

- `S00_raw`
- `S01_legacy_group2`
- `S02_run_only`
- `S03_minimal_run_plate`
- `S04_sequential_complete`
- `S05_sequential_complete_group_median`
- `S06_group_help`
- `S07_group_help_group_median`
- `S08_available_help`
- `S09_dual_anchor`
- `S10_batch_model`
- `S11_batch_help_pc_model`
- `S12_conditional_model`

脚本读取 `data/` 中的原始文件，按照仓库 README 中的名单剔除研究样本，并删除 Blank 与 Neat 控制。输出统一写入 `Results/`。

主入口 `run_all.R` 会依次执行 `script/parts/` 下的模块：

1. 函数和参数；
2. 数据载入与审计；
3. HELP 筛选及因子估计；
4. 多种 PACS 策略；
5. 统一评价；
6. 图形和自动报告。
