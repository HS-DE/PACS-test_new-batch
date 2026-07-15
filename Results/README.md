# PACS benchmark results

本目录由 `script/run_all.R` 自动生成。

输出结构：

- `audit/`：样本匹配、检出率、HELP 质量和实验设计检查；
- `factors/`：Run、Plate、样本级 HELP 因子及 HELP 主成分；
- `metrics/`：跨策略汇总和逐蛋白技术效应指标；
- `plots/`：PCA、校正因子及策略比较图；
- `matrices/`：各策略校正后的研究样本蛋白矩阵，压缩 RDS 格式；
- `REPORT.md`：自动生成的结果摘要；
- `strategy_manifest.tsv`：每种策略的准确定义。

生成文件不建议手工修改。需要修改方法时，应更新 `script/` 后重新运行分析。
