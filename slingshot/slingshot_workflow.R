# Slingshot trajectory workflow: edit the parameter block, then run with Rscript.
suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(slingshot)
})

input_rds <- "path/to/seurat_object.rds"
output_dir <- "results/04_trajectory/slingshot"
object_dir <- "data/processed"
assay <- "RNA"
reduction <- "umap"
cluster_col <- "cell_type"
start_cluster <- "RootCellType"
seed <- 123

set.seed(seed)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)

seu <- readRDS(input_rds)
DefaultAssay(seu) <- assay
counts <- GetAssayData(seu, assay = assay, layer = "counts")
metadata <- seu@meta.data[colnames(counts), , drop = FALSE]
embedding <- Embeddings(seu, reduction = reduction)[colnames(counts), 1:2, drop = FALSE]
metadata[[cluster_col]] <- factor(metadata[[cluster_col]])
if (!(start_cluster %in% levels(metadata[[cluster_col]]))) {
  stop("start_cluster not found in ", cluster_col, ": ", start_cluster, call. = FALSE)
}

sce <- SingleCellExperiment(
  assays = list(counts = counts),
  colData = metadata
)
reducedDims(sce)[[toupper(reduction)]] <- embedding
sce <- slingshot(
  sce,
  clusterLabels = cluster_col,
  reducedDim = toupper(reduction),
  start.clus = start_cluster
)

pseudotime <- slingPseudotime(sce, na = FALSE)
weights <- slingCurveWeights(sce)
write.csv(
  data.frame(cell = colnames(sce), pseudotime, check.names = FALSE),
  file.path(output_dir, "slingshot_pseudotime.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(cell = colnames(sce), weights, check.names = FALSE),
  file.path(output_dir, "slingshot_lineage_weights.csv"),
  row.names = FALSE
)

cluster_factor <- colData(sce)[[cluster_col]]
cluster_colors <- setNames(hcl.colors(nlevels(cluster_factor), "Dark 3"), levels(cluster_factor))
plot_trajectory <- function() {
  plot(
    embedding,
    col = cluster_colors[cluster_factor],
    pch = 16,
    cex = 0.5,
    xlab = colnames(embedding)[1],
    ylab = colnames(embedding)[2]
  )
  lines(SlingshotDataSet(sce), lwd = 2, col = "black")
  legend("topright", legend = levels(cluster_factor), col = cluster_colors, pch = 16, cex = 0.7, bty = "n")
}
pdf(file.path(output_dir, "slingshot_trajectory.pdf"), width = 8, height = 7)
plot_trajectory()
dev.off()
png(file.path(output_dir, "slingshot_trajectory.png"), width = 8, height = 7, units = "in", res = 300)
plot_trajectory()
dev.off()

object_path <- file.path(object_dir, "slingshot_sce.rds")
saveRDS(sce, object_path)
writeLines(c(
  "# Slingshot 轨迹分析结果",
  "",
  paste0("- 输入对象：", input_rds),
  paste0("- 聚类列：", cluster_col),
  paste0("- 起始 cluster：", start_cluster),
  paste0("- lineage 数量：", ncol(pseudotime)),
  paste0("- Slingshot SCE：", object_path)
), file.path(output_dir, "README_results_CN.md"))
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))