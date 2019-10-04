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

# old code for graph structure

cor_data <- as.data.frame(t(mydata))

# calculate correlation coefficients between all nodes
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
  mutate(adj.p.value = p.adjust(p.value, method = "BH")) %>%
  select(dim_x, everything())

square_cor_test <- lapply(cor_data, function(som_code){
  som_code_corr <- lapply(cor_data, cor.test, y = som_code) %>%
    lapply(broom::tidy) %>%
    bind_rows(.id = "dim_y")
  return(som_code_corr)
})

# # create a correlation distance matrix
# graph <- cor_test_results %>% 
#   mutate(weight = 1 - estimate,
#     V1 = dim_x,
#     V2 = dim_y) %>% 
#   select(V1,V2,weight) %>% 
#   graph_from_data_frame(directed = F)

# mst using euclidean distance
mst_euc <- mydata %>% 
  stats::dist(method = "euclidean") %>% 
  as.matrix() %>% 
  igraph::graph.adjacency(mode = "undirected",
                          weighted = T) %>% 
  igraph::minimum.spanning.tree()

mst_euc.layout <- igraph::layout_with_kk(mst_euc)

mst_euc.communities <- edge.betweenness.community(mst_euc,
                                                  weights=NULL,
                                                  directed=FALSE)
mst_euc.clustering <- make_clusters(mst_euc,
                                    membership=mst_euc.communities$membership)
V(mst_euc)$color <- mst_euc.communities$membership + 1

plot(mst_euc.clustering,
     mst_euc,
     layout = mst_euc.layout,
     vertex.label = NA)

# mst using pearson distance
mst_pearson <- mydata %>% 
  cor_dist("pearson") %>% 
  igraph::graph.adjacency(mode = "undirected",
                          weighted = T) %>% 
  igraph::minimum.spanning.tree()

mst_pearson_layout <- igraph::layout_with_kk(mst_pearson)

mst_pearson.communities <- edge.betweenness.community(mst_pearson, weights=NULL, directed=FALSE)
mst_pearson.clustering <- make_clusters(mst_pearson, membership=mst_pearson.communities$membership)
V(mst_pearson)$color <- mst_pearson.communities$membership + 1

plot(mst_pearson.clustering, mst_pearson, layout = mst_pearson_layout)

cluster_info <- data.frame(Barcode = rownames(reducedDim(sce_glm_pca, "GLM_PCA")),
                           som_id = som_model$unit.classif) %>% 
  left_join(data.frame(som_id = c(1:length(mst_pearson.communities$membership)),
                       cluster_id = mst_pearson.communities$membership))

# load in signatures

ar_reference_cpm <- read_csv("~/lgcp/rnaseq/analysis-scripts/ar_reference_beltran_cpm_grch37.csv")

neuro_reference_cpm <- read_csv("~/lgcp/rnaseq/analysis-scripts/neuro_reference_beltran_cpm_grch37.csv")

neuro <- read_csv("~/lgcp/rnaseq/analysis-scripts/neuro_reference_vpca_loadings_grch37.csv")

ar_signature_weights <- read_csv("~/lgcp/rnaseq/analysis-scripts/ar_signature_weights_Mendiratta.csv")

logcounts(sce_glm_pca)[rownames(logcounts(sce_glm_pca)) %in% ar_reference_cpm$Geneid,] %>% 
  as.matrix() %>% 
  hist()

counts(sce_glm_pca)[rownames(logcounts(sce_glm_pca)) %in% neuro_reference_cpm$Geneid,] %>% 
  as.matrix() %>% 
  hist()

logcounts(sce_glm_pca)[rowData(sce_glm_pca)$Symbol %in% neuro$Loading,] %>% 
  as.matrix() %>% 
  heatmap()

logcounts(sce_glm_pca)[rownames(logcounts(sce_glm_pca)) %in% ar_signature_weights$`Ensemble ID`,] %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  lapply(cor, y = ar_signature_weights$`BinReg Coef`)

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


barcode_ar_signature_score_count <- vector(length = ncol(sce_glm_pca))
for(i in 1:ncol(sce_glm_pca)){
  scored <- data.frame(logcounts = counts(sce_glm_pca)[rownames(counts(sce_glm_pca)) %in% ar_signature_weights$`Ensemble ID`,i],
                       ensgene = names(counts(sce_glm_pca)[rownames(counts(sce_glm_pca)) %in% ar_signature_weights$`Ensemble ID`,i])) %>% 
    left_join(ar_signature_weights, by = c("ensgene" = "Ensemble ID")) %>% 
    transmute(score = logcounts*`BinReg Coef`)
  
  barcode_ar_signature_score_count[i] <- sum(scored$score)
}

barcode_neuro_score_count <- vector(length = ncol(sce_glm_pca))
for(i in 1:ncol(sce_glm_pca)){
  scored <- data.frame(logcounts = counts(sce_glm_pca)[rowData(sce_glm_pca)$Symbol %in% neuro$Loading,i],
                       ensgene = names(counts(sce_glm_pca)[rowData(sce_glm_pca)$Symbol %in% neuro$Loading,i])) %>% 
    left_join(neuro_joined, by = c("ensgene" = "ID")) %>% 
    transmute(score = logcounts*PC1.v)
  
  barcode_neuro_score_count[i] <- sum(scored$score)
}

colData(sce_glm_pca)$som_id <- cluster_info$som_id
colData(sce_glm_pca)$cluster_id <- cluster_info$cluster_id

dgelist <- counts(sce_glm_pca) %>% 
  broom::tidy() %>% 
  setNames(c("ensgene", "Barcode", "count")) %>% 
  left_join(cluster_info) %>%  
  group_by(ensgene,cluster_id) %>% 
  summarise(count = sum(count)) %>% 
  spread(cluster_id, count, fill = 0) %>% 
  column_to_rownames("ensgene") %>% 
  as.matrix() %>% 
  edgeR::DGEList() %>% 
  calcNormFactors()

scored_clusters <- cpm(dgelist, log = T) %>% 
  as.data.frame() %>% 
  rownames_to_column("ensgene") %>% 
  left_join(ar_signature_weights %>% 
              transmute(ar_weight = `BinReg Coef`,
                        ensgene = `Ensemble ID`)) %>% 
  left_join(neuro_joined %>% 
              transmute(neuro_weight = PC1.v,
                        ensgene = ID)) %>% 
  left_join(neuro_reference_cpm %>% 
              transmute(neuro_reference = log2_cpm_mean,
                        ensgene = Geneid)) %>% 
  left_join(ar_reference_cpm %>% 
              setNames(c("ensgene", "ar_reference"))) %>% 
  gather(cluster_id,
         cpm,
         `1`:`30`) %>% 
  group_by(cluster_id) %>% 
  summarise(ar_signature_score = sum(ar_weight*cpm, na.rm = T),
            neuro_signature_score = sum(neuro_weight*cpm, na.rm = T),
            ar_correlation = cor(ar_reference, cpm, use = "pairwise.complete.obs", method = "spearman"),
            neuro_correlation = cor(neuro_reference, cpm, use = "pairwise.complete.obs", method = "spearman"))

mst_pearson_layout %>% 
  as.data.frame() %>% 
  setNames(c("mst_x", "mst_y")) %>% 
  rownames_to_column("som_id") %>%
  left_join(cluster_info %>% 
              distinct(som_id, cluster_id) %>% 
              mutate(som_id = as.character(som_id),
                     cluster_id = as.character(cluster_id))) %>% 
  left_join(scored_clusters) %>% 
  ggplot(aes(mst_x, mst_y, color = neuro_signature_score)) +
  geom_point() +
  theme_void()


save.image("/Volumes/group05/CCBB/CS024892_Kelly_Beshiri/ct-35-v1.RData")


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


plotReducedDim(sce, use_dimred = "UMAP")
plotReducedDim(sce, use_dimred = "PCA")

plotReducedDim(sce_glm_pca, use_dimred = "UMAP")
plotReducedDim(sce_glm_pca, use_dimred = "GLM_PCA", ncomponents = 5)
plotReducedDim(sce_glm_pca, use_dimred = "GLM_PCA", ncomponents = 6:10)