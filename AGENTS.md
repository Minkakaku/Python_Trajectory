# 轨迹分析项目约定

- 本仓库是独立分析点；六个根目录脚本分别接收 RDS/H5AD。
- root/start cluster 必须由生物学先验显式指定，不自动选择最方便的方向。
- PAGA 无 root 时只解释连通拓扑；scVelo 必须有真实 spliced/unspliced layers。
- tradeSeq 检验表保存 BH qvalue。
- 图保存 PNG/PDF，表格和 README_results_CN.md 放同一输出目录。
- legacy/original/ 和 source_catalog/ 保持原字节。