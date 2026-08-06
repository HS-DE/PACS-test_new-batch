# Adaptive HELP normalization v2

运行：

```bash
python -m pip install -r script/adaptive_help_v2/requirements.txt
python script/adaptive_help_v2/run_adaptive_help_v2.py
```

输出在 `Results/adaptive_help_v2/`。

## 流程概览

1. 去除已确认异常样本，蛋白强度做 log2。
2. 不看 HC/T2DM 标签，自动筛选稳定且不重复的 Core HELP。
3. Internal-QC HELP 用于判断和校正 Run 偏移。
4. pooled-QC HELP 用于板级 detection 校正；pooled-QC endogenous proteins 用于板级前处理校正。
5. 研究样本根据 HELP 偏移强度、HELP 一致性及 HELP/蛋白方向一致性，自动决定不校正、25% 校正或 50% 校正。
6. 比较 7 个候选流程，以 QC 稳定、Plate/Run 偏差降低和生物变异保留来选择最终流程。
7. HC/T2DM 标签只用于最后验证组内聚集、组间分离和分类表现，不参与 normalization factor 计算。
