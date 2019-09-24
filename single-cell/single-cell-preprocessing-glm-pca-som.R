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
# highly variable gene approach (standard)

source("single-cell/single-cell-scater-preprocessing-functions.R")

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

reducedDim(sce, "GLM_PCA") <- as.matrix(glmpca_poi_30$factors)

sce_glm_pca <- runUMAP(sce, use_dimred = "GLM_PCA", pca = 30)

som_grid <- somgrid(xdim = 30, ydim=30, topo="hexagonal")
# Finally, train the SOM, options for the number of iterations,
# the learning rates, and the neighbourhood are available
som_model <- som(reducedDim(sce, "GLM_PCA"), 
                 grid=som_grid, 
                 rlen=10000, 
                 alpha=c(0.05,0.01), 
                 keep.data = TRUE )

mydata <- getCodes(som_model)

cor_data <- as.data.frame(t(mydata))

cor_test_results <- lapply(seq_along(1:(ncol(cor_data)-1)),
       function(col_index){
         dim_x <- cor_data[,col_index]
         cols_to_test_start <- col_index + 1
         cols_to_test_end <- ncol(cor_data)
         colnames_index <- cols_to_test_start:cols_to_test_end
         dims_y <- cor_data[,cols_to_test_start:cols_to_test_end]
         dims_y <- as.data.frame(dims_y)
         colnames(dims_y) <- colnames(cor_data)[colnames_index]
         subs_cor_test_results <- lapply(dims_y, cor.test, y = dim_x) %>% 
           lapply(broom::tidy) %>% 
           bind_rows(.id = "dim_y") %>% 
           mutate(dim_x = rep(colnames(cor_data)[col_index], nrow(.)))
         return(subs_cor_test_results)
       }) %>% 
  bind_rows() %>% 
  mutate(adj.p.value = p.adjust(p.value)) %>% 
  select(dim_x, everything())

# create a correlation matrix where insignificant correlations are set to 0, then do the graphing
graph <- cor_test_results %>% 
  mutate(weight = if_else(estimate > 0.5, estimate, 0),
         V1 = dim_x,
         V2 = dim_y) %>% 
  select(V1,V2,weight) %>% 
  graph_from_data_frame(directed = F)

plot.igraph(graph, edge.width = E(graph)$edge.width, 
            edge.color = "blue", vertex.color = "white", vertex.size = 1,
            vertex.frame.color = NA, vertex.label.color = "grey30")
clp <- cluster_walktrap(graph)

som_model$unit.classif

clp$membership

cluster_info <- data.frame(Barcode = rownames(reducedDim(sce_glm_pca, "GLM_PCA")),
                           som_id = som_model$unit.classif) %>% 
  left_join(data.frame(som_id = c(1:length(clp$membership)),
                       cluster_id = clp$membership))

deg <- findMarkers(sce_glm_pca,
                   cluster_info$cluster_id,
                   direction = "up",
                   pval.type="all") %>% 
  lapply(as.data.frame) %>% 
  lapply(rownames_to_column, "ensgene") %>% 
  bind_rows(.id = "cluster_id") %>% 
  left_join(grch38)

deg %>% 
  dplyr::filter(FDR < 0.05) %>% 
  View()

sce_glm_pca$Clusters = factor(cluster_info$cluster_id)

save.image("/Volumes/group05/CCBB/CS024892_Kelly_Beshiri/Untitled-ct-35.RData")
