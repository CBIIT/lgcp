if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

BiocManager::install(c("DropletUtils"))


library(DropletUtils)
library(scater)

sce <- read10xCounts("/Volumes/Group05/CCBB/SA/Bioinfo_2018/CS024650_Kelly_Agarwal/02_PrimaryAnalysisOutput/00_FullCellrangerOutput/SCAF604_683/outs/filtered_feature_bc_matrix/")
ctrl <- list(Mito = grep("^MT|^mt-", rowData(sce)$Symbol)) 

sce <- calculateQCMetrics(sce,
                          feature_controls = ctrl,
                          compact = TRUE,
                          BPPARAM = MulticoreParam(7))

high_mito <- isOutlier(colData(sce)$scater_qc$feature_control_Mito$pct_counts,
                       nmads = 3, type = "higher")
metadata(sce)$high_mito <- high_mito
sce <- sce[, !metadata(sce)$high_mito]

num_reads <- 1
num_cells <- 0.05 * ncol(sce)
keep <- which(DelayedArray::rowSums(counts(sce) >= num_reads) >= num_cells)

metadata(sce)$genes_keep <- keep

low_lib_sce <- isOutlier(sce$log10_total_counts, type="lower", nmad=3)
low_genes_sce <- isOutlier(sce$log10_total_features_by_counts, type="lower", nmad=3)

sce <- sce[, !(low_lib_sce | low_genes_sce)]
