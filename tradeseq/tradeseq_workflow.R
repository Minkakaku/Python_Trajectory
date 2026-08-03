# tradeSeq downstream analysis for a Slingshot SingleCellExperiment.
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(slingshot)
  library(tradeSeq)
})

input_rds <- "data/processed/slingshot_sce.rds"
output_dir <- "results/04_trajectory/tradeseq"
object_dir <- "data/processed"
nknots <- 6
workers <- 4

bpparam <- if (.Platform$OS.type == "windows") {
  BiocParallel::SnowParam(workers = workers)
} else {
  BiocParallel::MulticoreParam(workers = workers)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)

sce <- readRDS(input_rds)
counts <- assay(sce, "counts")
pseudotime <- slingPseudotime(sce, na = FALSE)
cell_weights <- slingCurveWeights(sce)

set.seed(123)
gam <- fitGAM(
  counts = counts,
  pseudotime = pseudotime,
  cellWeights = cell_weights,
  nknots = nknots,
  parallel = workers > 1,
  BPPARAM = bpparam,
  verbose = TRUE
)

association <- associationTest(gam)
association$gene <- rownames(association)
association <- association[order(association$pvalue), , drop = FALSE]
write.csv(association, file.path(output_dir, "association_test.csv"), row.names = FALSE)

start_end <- startVsEndTest(gam)
start_end$gene <- rownames(start_end)
start_end <- start_end[order(start_end$pvalue), , drop = FALSE]
write.csv(start_end, file.path(output_dir, "start_vs_end_test.csv"), row.names = FALSE)

if (ncol(cell_weights) > 1) {
  pattern <- patternTest(gam)
  pattern$gene <- rownames(pattern)
  pattern <- pattern[order(pattern$pvalue), , drop = FALSE]
  write.csv(pattern, file.path(output_dir, "pattern_test.csv"), row.names = FALSE)

  diff_end <- diffEndTest(gam)
  diff_end$gene <- rownames(diff_end)
  diff_end <- diff_end[order(diff_end$pvalue), , drop = FALSE]
  write.csv(diff_end, file.path(output_dir, "differential_end_test.csv"), row.names = FALSE)
}

object_path <- file.path(object_dir, "tradeseq_gam.rds")
saveRDS(gam, object_path)
writeLines(c(
  "# tradeSeq 轨迹差异表达结果",
  "",
  paste0("- 输入 Slingshot 对象：", input_rds),
  paste0("- lineage 数量：", ncol(cell_weights)),
  paste0("- nknots：", nknots),
  paste0("- tradeSeq GAM：", object_path)
), file.path(output_dir, "README_results_CN.md"))
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))