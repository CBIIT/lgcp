library(DropletUtils)
library(scater)
library(scran)
library(annotables)
library(tidyverse)
library(flowCore)
library(FlowSOM)
library(glmpca)

# highly variable gene approach (standard)

sce <- read10xCounts("/Volumes/group09/CCB/Beshiri/Folders_old/CT35/'Omics_data/single_cell_RNAseq/CS024464_Beshiri_CellTagging/CS024464_Beshiri_CellTagging/02-CountsOutput/SCAF636_35-1/outs/filtered_feature_bc_matrix/")

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
hvg <- rownames(dec[dec$bio > 0, ]) # ~4k genes

## dimensionality reduction
set.seed(1234)
sce <- runPCA(sce, 
              feature_set = hvg)

sce <- runUMAP(sce)

# setwd("/Volumes/Home04/CapaldoBj/scrna2019-master/")
# source("/Volumes/Home04/CapaldoBj/scrna2019-master/util/functions_genefilter.R")

colnames(sce) <- colData(sce)$Barcode
ranked_genes <- rank_all_genes(sce, "total_counts")

sce_d <- sce[ranked_genes$dev <= 2000,]

filtered_counts <- counts(sce_d)
filtered_counts <- filtered_counts[rowSums(filtered_counts) > 0,]

glmpca_poi_30 <- glmpca(as.matrix(filtered_counts),
                             30,
                             fam = "poi")
save.image("/Volumes/group05/CCBB/CS024892_Kelly_Beshiri/Untitled-ct-35.RData")

reducedDim(sce, "GLM_PCA") <- as.matrix(glmpca_poi_30$factors)

sce_glm_pca <- runUMAP(sce, use_dimred = "GLM_PCA", pca = 30)

g <- buildSNNGraph(sce_glm_pca, k=25, use.dimred = 'GLM_PCA', type="rank")
clust <- igraph::cluster_louvain(g)$membership
table(clust)

plotReducedDim(sce_glm_pca, "UMAP", colour_by="cluster")
plotReducedDim(sce_glm_pca, "GLM_PCA", colour_by="cluster")
plotReducedDim(sce_glm_pca, "GLM_PCA", ncomponents = 5, colour_by="cluster")
plotReducedDim(sce_glm_pca, "GLM_PCA", ncomponents = 6:10, colour_by="cluster")


flow_frame <- new("flowFrame",
                  exprs = as.matrix(glmpca_poi_30$factors))

som <- ReadInput(flow_frame)
som <- BuildSOM(som, xdim = 5, ydim = 5)
mst <- BuildMST(som)

cluster_diff <- scran::findMarkers(sce, mst$map$mapping[,1])

cluster_diff %>% 
  lapply(as.data.frame) %>% 
  bind_rows(.id = "cluster_id") %>% 
  mutate(cluster_id = factor(cluster_id, levels = c(1:25))) %>% 
  ggplot(aes(p.value)) + 
  geom_histogram(binwidth = 0.05) +
  facet_wrap(~cluster_id) +
  theme_bw()