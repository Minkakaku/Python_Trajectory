# R-Python-trajectory-workflow-template：轨迹、拟时序与 RNA velocity

这里集中保存所有能用于单细胞轨迹分析的代码。新项目优先从六个正式入口选择方法；历史项目、平台封装及其他人员代码保留在来源层，用于查阅实现细节和补充特殊分析。

## 目录结构

```text
R-Python-trajectory-workflow-template/
├── monocle2/                  # 正式 Monocle2 + BEAM 入口
├── monocle3/                  # 正式 Monocle3 principal graph 入口
├── slingshot/                 # 正式 Slingshot lineage 入口
├── tradeseq/                  # 正式 tradeSeq 轨迹差异表达入口
├── paga/                      # 正式 Scanpy PAGA 入口
├── rna_velocity/              # 正式 scVelo RNA velocity 入口
├── source_catalog/            # 343 份轨迹相关代码的逐字节副本
├── source_archive/            # 原项目/平台脚本，保留硬编码与原始写法
└── SOURCE_CATALOG.tsv         # 方法、场景、来源路径、哈希和目标文件索引
```

- 正式入口已修正明显错误，路径可配置，图表和表格写入 `results/04_trajectory/<method>/`。
- 中间对象写入 `data/processed/`，不放进 `results/`。
- `source_catalog/` 与 `../R-Python-bioinformatics-code-archive/11_semantic_catalog/` 中对应文件 SHA-256 完全一致，不改原代码。
- `source_archive/` 保留原平台、Streamlit、Docker、项目路径与样本设置，只用于溯源，不建议直接作为新项目入口。

## 怎么选方法

| 方法 | 适合的问题 | 主要输入 | 说明 |
|---|---|---|---|
| Monocle2 | 经典 DDRTree 拟时序、分支与 BEAM | Seurat raw counts、聚类/注释列 | 兼容旧项目多，BEAM 分支分析完整；依赖较旧 |
| Monocle3 | principal graph、partition 内轨迹和 graph_test | Seurat raw counts、UMAP、起始细胞类型 | 图结构更灵活，可指定 root cell group |
| Slingshot | 已有聚类和降维上的 lineage 推断 | Seurat、UMAP、起始 cluster | 不重新定义表达空间，适合稳定聚类后的轨迹 |
| tradeSeq | 沿 Slingshot lineage 的动态差异表达 | `slingshot_sce.rds` | 提供 association、start-vs-end、pattern、diffEnd 检验 |
| PAGA | cluster 间连通性和粗粒度拓扑 | AnnData、分类分组列、邻接图 | 没有邻接图时正式脚本会自动计算 neighbors |
| scVelo | 方向性 RNA velocity 与 latent time | 含 `spliced`/`unspliced` layers 的 h5ad | 支持 deterministic、stochastic、dynamical 模式 |
| CytoTRACE 等 | 分化潜能或其他轨迹辅助判断 | 依原代码而定 | 目前保留在 `source_catalog/`，按索引挑选 |

Monocle/PAGA/Slingshot主要推断状态连续性或拓扑；scVelo利用剪接动力学提供方向信息；tradeSeq是在已经得到的 lineage 上检验基因动态。它们可以互相验证，但结果含义并不相同。

## 正式入口

### Monocle2

```bash
Rscript monocle2/monocle2_workflow.R \
  --rds data/processed/seurat.rds \
  --output results/04_trajectory/monocle2 \
  --object_dir data/processed \
  --cluster cell_type \
  --root_col cell_type \
  --root_group RootCellType \
  --BEAM
```

`--downcluster --maxcells 1000` 可按 cluster 下采样。正式脚本保留了原代码的 dpFeature、DDRTree、root state 与 BEAM 逻辑，并修正普通矩阵检查、branch point 数量及 BEAM 输出缺少 gene 列的问题。

### Monocle3

编辑 `monocle3/monocle3_workflow.R` 顶部参数后运行：

```bash
Rscript monocle3/monocle3_workflow.R
```

把 `root_col` 和 `root_groups` 设置为生物学上合理的起始细胞。脚本会把这些细胞投票到 principal graph 上最常见的节点作为 root。

### Slingshot 与 tradeSeq

先编辑并运行 Slingshot，再运行 tradeSeq：

```bash
Rscript slingshot/slingshot_workflow.R
Rscript tradeseq/tradeseq_workflow.R
```

Slingshot 的 `start_cluster` 必须是 `cluster_col` 中存在的标签。tradeSeq 默认读取 `data/processed/slingshot_sce.rds`，并在 Windows 使用 SnowParam、Linux/macOS 使用 MulticoreParam。

### PAGA

```bash
python paga/paga_workflow.py \
  --input-h5ad data/processed/adata.h5ad \
  --groupby cell_type
```

可用 `--threshold`、`--n-pcs`、`--n-neighbors` 调整图连通阈值和邻接图参数。

### scVelo

```bash
python rna_velocity/scvelo_workflow.py \
  --input-h5ad data/processed/velocity.h5ad \
  --groupby cell_type \
  --basis umap \
  --mode dynamical
```

dynamical 模式会额外运行 `recover_dynamics()` 并输出 latent time；数据必须包含 spliced 与 unspliced layers。

## 输出规范

每个正式入口使用自己的结果目录：

```text
results/04_trajectory/<method>/
├── *.pdf
├── *.png
├── *.csv 或 *.xls
├── README_results_CN.md
└── sessionInfo.txt            # R 入口

data/processed/
├── monocle2_cds.rds
├── monocle3_cds.rds
├── slingshot_sce.rds
├── tradeseq_gam.rds
├── paga.h5ad
└── scvelo.h5ad
```

正式图均保存同名 PDF 和 300-dpi PNG；伪时间、lineage weights、连通矩阵和轨迹检验结果保留为表格。

## 全量来源索引

`SOURCE_CATALOG.tsv` 当前记录 343 份唯一 SHA-256 代码，来源包括：

- 281 份 Monocle 专题代码；
- 21 份其他轨迹专题代码；
- 31 份嵌在核心单细胞流程中的轨迹代码；
- 10 份被聚类、富集或 Xenium 分类收纳、但实际调用轨迹方法的代码。

同一文件可同时带多个 `method` 标签，例如 `monocle2;monocle2_BEAM`。查找特殊实现时，优先按 `method`、`analysis_scenario`、`workflow_role` 筛选，不要按原文件名去重；同名脚本可能分别负责主流程、分支分析、绘图或平台调用。

## 使用原则

1. 新项目从正式入口或 `source_catalog/` 选择代码，组合成新的 `workflow/stepXX.分析点.R` 或 `.py`。
2. 不直接修改已经运行过的 workflow step；重跑时新增下一编号。
3. root、细胞类型列、样本分组和参数属于项目决策，应在新项目配置或 workflow 中明确记录。
4. 自动轨迹不能替代生物学方向判断；root 选择、已知 marker、采样时间和不同方法的一致性需要一起检查。
5. 本次完成的是代码整理与静态语法检查，没有用真实表达对象完整复跑六种方法。
