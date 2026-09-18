# R-Python-trajectory-workflow-template

独立轨迹分析仓库。六种方法各保留一个根目录完整脚本：

- monocle2_workflow.R：DDRTree、dpFeature ordering genes、可选 BEAM。
- monocle3_workflow.R：principal graph、显式 root、graph_test。
- slingshot_workflow.R：cluster-level lineage 与 principal curves。
- tradeseq_workflow.R：可直接接收 Seurat RDS 或 Slingshot SCE，拟合 NB-GAM。
- paga_workflow.py：PAGA 群间拓扑；指定 root 后可计算 DPT。
- scvelo_workflow.py：基于 spliced/unspliced 的 RNA velocity。

## 方法选择

- Monocle2 用于复现旧项目；新分析优先考虑 Monocle3、Slingshot/tradeSeq 或 PAGA。
- Monocle3 和 Slingshot 必须显式指定有生物学依据的起始群。
- tradeSeq 用于检验沿 pseudotime 的表达变化及多 lineage 差异，不负责决定 root。
- PAGA 本身只表示群间连通结构；没有 root 时不能解释为有方向的发育顺序。
- scVelo 只有在 h5ad 保留真实 spliced/unspliced layers 时才能运行，普通表达矩阵不能补造 RNA velocity。

## 示例

    Rscript monocle3_workflow.R --input annotated.rds --output results/monocle3 --root-col cell_type --root-groups Progenitor

    Rscript slingshot_workflow.R --input annotated.rds --output results/slingshot --cluster-col cell_type --start-cluster Progenitor

    Rscript tradeseq_workflow.R --input annotated.rds --output results/tradeseq --input-type seurat --start-cluster Progenitor

    python paga_workflow.py --input-h5ad annotated.h5ad --output-dir results/paga --groupby cell_type --root-group Progenitor

    python scvelo_workflow.py --input-h5ad velocity.h5ad --output-dir results/scvelo --groupby cell_type --mode stochastic

所有 pseudotime/latent time 均为模型中的相对顺序，不等于真实时间。结果需结合采样设计、marker、已知起终点和其他方法复核。

legacy/original/ 保存来源脚本；source_catalog/ 保存语义副本，不作为正式入口。

方法依据：Monocle2（PMID: 28825705）、Monocle3 官方文档、Slingshot（PMID: 29914354）、tradeSeq（PMID: 32139671）、PAGA（PMID: 30890159）、scVelo（PMID: 32747759）。