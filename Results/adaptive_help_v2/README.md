# Adaptive HELP normalization v2 — 运行结果

## 最终选择

**S4_adaptive_sample_help**

本次流程没有使用 HC/T2DM 标签计算任何 normalization factor。标签只用于最后检查组内聚集和组间分离。

## 主要处理

1. 去除 README 中已确认的 29 个异常样本，并对蛋白强度做 log2。
2. 从 50 条 HELP 中，按检出率、Internal-QC 稳定性、pooled-QC 稳定性和冗余程度，自动选择 25 条 Core HELP。
3. 用 Internal-QC HELP 判断并校正仪器 Run 偏移。
4. 用 pooled-QC HELP 校正板级 detection 偏移，再用 pooled-QC endogenous proteins 校正板级前处理偏移。
5. 每个研究样本先在同一 Run/Plate 区块内判断 HELP 偏移是否真实且与整体蛋白偏移一致，再决定 0%、25% 或 50% 校正；不再对所有样本一刀切。
6. 比较 7 个候选流程，以 QC 变异、Plate/Run 残留、HELP/蛋白残留相关和生物变异保留作为主选择指标；HC/T2DM 聚集和 AUC 只作验证。

## 关键结果

- Study samples: **185**
- Core HELP: **25**
- Run correction reliability: **0.118**
- Plate HELP correction reliability: **0.569**
- 25% sample HELP correction: **9 samples**
- 50% sample HELP correction: **14 samples**
- QC median log2 SD: **0.0642**
- QC SD improvement: **14.6%**
- Plate eta²: **0.0400**
- Run eta²: **0.0089**
- Within-group distance: **17.4219**
- Within-group improvement vs raw: **0.16%**
- Between/within ratio: **0.4861**
- Between/within improvement vs raw: **1.87%**
- Repeated CV AUC: **0.602 ± 0.054**
- Leave-one-plate-out AUC: **0.682**

完整数字见 `metrics/strategy_summary.csv`，最终矩阵见 `matrices/final_normalized_study_matrix.tsv.gz`。
