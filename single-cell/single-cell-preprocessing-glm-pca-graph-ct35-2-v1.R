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

rowData(sce_glm_pca) <- cbind(rowData(sce_glm_pca),
                              ranked_genes)

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

# perform som-kmeans-signature-scoring

full_som_grid <- somgrid(xdim = 30, ydim=30, topo="rectangular")

full_som_model <- som(reducedDim(sce, "GLM_PCA"), 
                      grid=full_som_grid, 
                      rlen=10000, 
                      alpha=c(0.05,0.01), 
                      keep.data = TRUE )

som_grid <- somgrid(xdim = 2, ydim=2, topo="rectangular")
# Finally, train the SOM, options for the number of iterations,
# the learning rates, and the neighbourhood are available

som_results <- lapply(seq_along(1:(ncol(reducedDim(sce, "GLM_PCA"))-1)),
                      function(col_index){
                        dim_x <- reducedDim(sce, "GLM_PCA")[,col_index]
                        cols_to_test_start <- col_index + 1
                        cols_to_test_end <- ncol(reducedDim(sce, "GLM_PCA"))
                        colnames_index <- cols_to_test_start:cols_to_test_end
                        dims_y <- reducedDim(sce, "GLM_PCA")[,cols_to_test_start:cols_to_test_end]
                        dims_y <- as.data.frame(dims_y)
                        colnames(dims_y) <- colnames(reducedDim(sce, "GLM_PCA"))[colnames_index]
                        
                        subs_som_results <- lapply(dims_y,
                                                   function(dim_y){
                                                     som_data <- data.frame(dim_x = dim_x,
                                                                            dim_y = dim_y) %>% 
                                                       as.matrix()
                                                     som_model <- som(som_data, 
                                                                      grid=som_grid, 
                                                                      rlen=10000, 
                                                                      alpha=c(0.05,0.01), 
                                                                      keep.data = TRUE )
                                                     return(som_model$unit.classif)
                                                   })
                        return(subs_som_results)
                      }) %>% 
  lapply(bind_cols) %>%
  lapply(rownames_to_column, "cell_index")
names(som_results) <- paste0("som_dim", 1:length(som_results))

som_results <- som_results %>% 
  bind_rows(.id = "dim_x") %>% 
  gather(dim_y,
         som_code,
         -c(dim_x, cell_index)) %>% 
  filter(!is.na(som_code)) %>% 
  mutate(som_code = factor(som_code)) %>% 
  unite(dim_x_y,
        dim_x,
        dim_y) %>% 
  spread(dim_x_y,
         som_code)

kmeans_results <- lapply(seq_along(1:(ncol(reducedDim(sce, "GLM_PCA"))-1)),
                         function(col_index){
                           dim_x <- reducedDim(sce, "GLM_PCA")[,col_index]
                           cols_to_test_start <- col_index + 1
                           cols_to_test_end <- ncol(reducedDim(sce, "GLM_PCA"))
                           colnames_index <- cols_to_test_start:cols_to_test_end
                           dims_y <- reducedDim(sce, "GLM_PCA")[,cols_to_test_start:cols_to_test_end]
                           dims_y <- as.data.frame(dims_y)
                           colnames(dims_y) <- colnames(reducedDim(sce, "GLM_PCA"))[colnames_index]
                           
                           subs_kmeans_results <- lapply(dims_y,
                                                         function(dim_y){
                                                           kmeans_data <- data.frame(dim_x = dim_x,
                                                                                     dim_y = dim_y) %>% 
                                                             as.matrix()
                                                           kmeans_clust <- kmeans(kmeans_data, centers = 4, iter.max = 1000)
                                                           return(kmeans_clust$cluster)
                                                         })
                           return(subs_kmeans_results)
                         }) %>% 
  lapply(bind_cols) %>%
  lapply(rownames_to_column, "cell_index")
names(kmeans_results) <- paste0("kmeans_dim", 1:length(kmeans_results))

kmeans_results <- kmeans_results %>% 
  bind_rows(.id = "dim_x") %>% 
  gather(dim_y,
         som_code,
         -c(dim_x, cell_index)) %>% 
  filter(!is.na(som_code)) %>% 
  mutate(som_code = factor(som_code)) %>% 
  unite(dim_x_y,
        dim_x,
        dim_y) %>% 
  spread(dim_x_y,
         som_code)

split_results <- lapply(seq_along(1:(ncol(reducedDim(sce, "GLM_PCA"))-1)),
                        function(col_index){
                          dim_x <- reducedDim(sce, "GLM_PCA")[,col_index]
                          cols_to_test_start <- col_index + 1
                          cols_to_test_end <- ncol(reducedDim(sce, "GLM_PCA"))
                          colnames_index <- cols_to_test_start:cols_to_test_end
                          dims_y <- reducedDim(sce, "GLM_PCA")[,cols_to_test_start:cols_to_test_end]
                          dims_y <- as.data.frame(dims_y)
                          colnames(dims_y) <- colnames(reducedDim(sce, "GLM_PCA"))[colnames_index]
                          
                          subs_kmeans_results <- lapply(dims_y,
                                                        function(dim_y){
                                                          kmeans_data <- data.frame(dim_x = dim_x,
                                                                                    dim_y = dim_y) %>% 
                                                            rownames_to_column("barcode") %>% 
                                                            gather(dim,
                                                                   position,
                                                                   starts_with("dim")) %>% 
                                                            group_by(barcode) %>% 
                                                            summarize(distance = (sum(position^2))^(1/2))
                                                          kmeans_clust <- kmeans(kmeans_data$distance, centers = 2, iter.max = 1000)
                                                          if(kmeans_clust$centers[1,1] > kmeans_clust$centers[2,1]){
                                                            clusters <- data.frame(clust_id = kmeans_clust$cluster) %>% 
                                                              mutate(clust_id = if_else(clust_id == 1, 1, 0))
                                                          }else{
                                                            clusters <- data.frame(clust_id = kmeans_clust$cluster) %>% 
                                                              mutate(clust_id = if_else(clust_id == 2, 1, 0))
                                                          }
                                                          return(clusters$clust_id)
                                                        })
                          return(subs_kmeans_results)
                        }) %>% 
  lapply(bind_cols) %>%
  lapply(rownames_to_column, "cell_index")
names(split_results) <- paste0("split_dim", 1:length(split_results))
split_results <- split_results %>% 
  bind_rows(.id = "dim_x") %>% 
  gather(dim_y,
         som_code,
         -c(dim_x, cell_index)) %>% 
  filter(!is.na(som_code)) %>% 
  mutate(som_code = factor(som_code)) %>% 
  unite(dim_x_y,
        dim_x,
        dim_y) %>% 
  spread(dim_x_y,
         som_code)

colData(sce_glm_pca)$cell_index <- as.character(c(1:nrow(colData(sce_glm_pca))))
colData(sce_glm_pca)$full_som_node <- full_som_model$unit.classif

colData(sce_glm_pca) <- colData(sce_glm_pca) %>%
  as.data.frame() %>%
  left_join(som_results) %>%
  left_join(kmeans_results) %>%
  left_join(split_results) %>%
  DataFrame()

save(sce_glm_pca, file = "/Volumes/Group05/CCBB/Single-Cell-Bioinformatics-2019-October-03/ct-35-2-v1-sce-object.Rdata")

# NE score ranking
logcounts(sce_glm_pca)[rowData(sce_glm_pca)$Symbol %in% neuro$Loading,] %>% 
  tidy_sparse_matrix() %>% 
  left_join(colData(sce_glm_pca) %>% 
              as.data.frame() %>% 
              select(Barcode, clust),
            by = c("column" = "Barcode")) %>% 
  left_join(neuro_joined, by = c("row" = "ID")) %>% 
  mutate(weighted_expression = value*PC1.v) %>% 
  group_by(clust, Loading) %>% 
  summarize(mean_weighted_expression = mean(weighted_expression, na.rm = T)) %>% 
  spread(clust,
         mean_weighted_expression) %>% 
  write_csv("/Volumes/Group05/CCBB/Single-Cell-Bioinformatics-2019-October-03/NE-score-weighted-expression-values.csv")

# AR score ranking
logcounts(sce_glm_pca)[rownames(logcounts(sce_glm_pca)) %in% ar_signature_weights$`Ensemble ID`,] %>% 
  broom::tidy() %>% 
  left_join(colData(sce_glm_pca) %>% 
              as.data.frame() %>% 
              select(Barcode, clust),
            by = c("column" = "Barcode")) %>% 
  left_join(ar_signature_weights, by = c("row" = "Ensemble ID")) %>% 
  mutate(weighted_expression = value*`BinReg Coef`) %>% 
  group_by(clust, `Gene Symbol`) %>% 
  summarize(mean_weighted_expression = mean(weighted_expression, na.rm = T)) %>% 
  spread(clust,
         mean_weighted_expression) %>% 
  write_csv("/Volumes/Group05/CCBB/Single-Cell-Bioinformatics-2019-October-03/AR-signature-score-weighted-expression-values.csv")

reducedDim(sce_glm_pca, "UMAP") %>%
  as.data.frame() %>%
  rownames_to_column("Barcode") %>%
  setNames(c("Barcode", "UMAP-1", "UMAP-2")) %>%
  write_csv("/Volumes/Group05/CCBB/Single-Cell-Bioinformatics-2019-October-03/ct35-2-v1-UMAP-projection.csv")

colData(sce_glm_pca) %>%
  as.data.frame() %>%
  select(Barcode, clust) %>% 
  write_csv("/Volumes/Group05/CCBB/Single-Cell-Bioinformatics-2019-October-03/ct35-2-v1-glm-pca-graph-clusters.csv")

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

