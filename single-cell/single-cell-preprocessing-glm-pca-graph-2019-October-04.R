library(DropletUtils)
library(scater)
library(scran)
library(annotables)
library(tidyverse)
library(glmpca)
library(kohonen)
library(ggraph)
library(tidygraph)
library(igraph)
library(Seurat)
library(edgeR)

set.seed(8675309)

# highly variable gene approach (standard)

source("single-cell/single-cell-scater-preprocessing-functions.R")

# load in signatures
neuro <- read_csv("~/lgcp/rnaseq/analysis-scripts/neuro_reference_vpca_loadings_grch37.csv")
ar_signature_weights <- read_csv("~/lgcp/rnaseq/analysis-scripts/ar_signature_weights_Mendiratta.csv")


sce <- read10xCounts("/Volumes/Group09/CCB/Beshiri/Folders_old/CT35/'Omics_data/single_cell_RNAseq/lineage-tracing-May-2018/CT35_10x_filtered_gbm")

# ADD FEATURE FOR REF ORG HERE!
location_tidy <- rowData(sce) %>% 
  as.data.frame() %>% 
  left_join(grch38 %>% 
              select(ensgene, chr), by = c("ID" = "ensgene")) %>% 
  distinct()
is.mito <- which(location_tidy$chr == "MT")

sce <- calculateQCMetrics(sce, feature_controls=list(Mito=is.mito), BPPARAM = MulticoreParam(7))

# remove cells with log-library sizes that are more than 3 MADs below the median
qc.lib <- isOutlier(log(sce$total_counts), nmads=3, type="lower")
# remove cells where the log-transformed number of expressed genes is 3 MADs below the median
qc.nexprs <- isOutlier(log(sce$total_features_by_counts), nmads=3,
                       type="lower")
# remove cells where the number of mito genes is 3 MADs above the median
qc.mito <- isOutlier(sce$pct_counts_Mito, nmads=3, type="higher")

discard <- qc.lib | qc.nexprs | qc.mito

# Summarize the number of cells removed for each reason.
cells_removed_summary_df <- DataFrame(LibSize=sum(qc.lib),
                                      NExprs=sum(qc.nexprs),
                                      MitoProp=sum(qc.mito),
                                      Total=sum(discard))

# Retain only high-quality cells in the SingleCellExperiment.
sce <- sce[,!discard]

## normalization

cl<-scran::quickCluster(sce)
sce<-scran::computeSumFactors(sce,clusters=cl)
sce <- scater::normalize(sce)

## feature selection
fit <- trendVar(sce, use.spikes = F)
dec <- decomposeVar(sce, fit)

colnames(sce) <- colData(sce)$Barcode
ranked_genes <- rank_all_genes(sce, "total_counts")

sce_d <- sce[ranked_genes$dev <= 2000,]

filtered_counts <- counts(sce_d)
filtered_counts <- filtered_counts[rowSums(filtered_counts) > 0,]

glmpca_poi_30 <- glmpca(as.matrix(filtered_counts),
                        30,
                        fam = "poi")

reducedDim(sce, "GLM_PCA") <- as.matrix(glmpca_poi_30$factors)

sce_glm_pca <- runUMAP(sce, use_dimred = "GLM_PCA", pca = 30)

barcode_ar_signature_score <- vector(length = ncol(sce_glm_pca))
for(i in 1:ncol(sce_glm_pca)){
  scored <- data.frame(logcounts = logcounts(sce_glm_pca)[rownames(logcounts(sce_glm_pca)) %in% ar_signature_weights$`Ensemble ID`,i],
             ensgene = names(logcounts(sce_glm_pca)[rownames(logcounts(sce_glm_pca)) %in% ar_signature_weights$`Ensemble ID`,i])) %>% 
    left_join(ar_signature_weights, by = c("ensgene" = "Ensemble ID")) %>% 
    transmute(score = logcounts*`BinReg Coef`)
  
  barcode_ar_signature_score[i] <- sum(scored$score)
}

neuro_joined <- neuro %>% 
  left_join(as.data.frame(rowData(sce_glm_pca)), by = c("Loading" = "Symbol"))

barcode_neuro_score <- vector(length = ncol(sce_glm_pca))
for(i in 1:ncol(sce_glm_pca)){
  scored <- data.frame(logcounts = logcounts(sce_glm_pca)[rowData(sce_glm_pca)$Symbol %in% neuro$Loading,i],
                       ensgene = names(logcounts(sce_glm_pca)[rowData(sce_glm_pca)$Symbol %in% neuro$Loading,i])) %>% 
    left_join(neuro_joined, by = c("ensgene" = "ID")) %>% 
    transmute(score = logcounts*PC1.v)
  
  barcode_neuro_score[i] <- sum(scored$score)
}

# perform graph based clustering

g <- buildSNNGraph(sce_glm_pca, k=10, use.dimred = 'GLM_PCA')
clust <- igraph::cluster_walktrap(g)$membership
table(clust)

colData(sce_glm_pca)$clust <- factor(clust)
colData(sce_glm_pca)$NE_score <- barcode_neuro_score
colData(sce_glm_pca)$AR_score <- barcode_ar_signature_score

plotReducedDim(sce_glm_pca, use_dimred = "UMAP", colour_by = "clust")
plotReducedDim(sce_glm_pca, use_dimred = "UMAP", colour_by = "NE_score")
plotReducedDim(sce_glm_pca, use_dimred = "UMAP", colour_by = "AR_score")

deg_all <- findMarkers(sce_glm_pca,
                   sce_glm_pca$clust,
                   direction = "up",
                   pval.type="all") %>% 
  lapply(as.data.frame) %>% 
  lapply(rownames_to_column, "ensgene") %>% 
  bind_rows(.id = "cluster_id") %>% 
  left_join(grch38)

deg_any <- findMarkers(sce_glm_pca,
                       sce_glm_pca$clust,
                       direction = "up",
                       pval.type="any") %>% 
  lapply(as.data.frame) %>% 
  lapply(rownames_to_column, "ensgene") %>% 
  bind_rows(.id = "cluster_id") %>% 
  left_join(grch38)

for(clust_id in unique(deg_all$cluster_id)){
  file <- paste0("/Volumes/Group05/CCBB/Single-Cell-Bioinformatics-2019-October-03/", "ct-35-2-v1-10x-deg-all-cluster-id-", clust_id, ".csv")
  deg_all %>% 
    filter(cluster_id == clust_id) %>% 
    write_csv(file)
}

for(clust_id in unique(deg_any$cluster_id)){
  file <- paste0("/Volumes/Group05/CCBB/Single-Cell-Bioinformatics-2019-October-03/", "ct-35-2-v1-10x-deg-any-cluster-id-", clust_id, ".csv")
  deg_any %>% 
    filter(cluster_id == clust_id) %>% 
    write_csv(file)
}

