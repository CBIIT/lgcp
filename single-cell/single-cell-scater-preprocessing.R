library(DropletUtils)
library(scater)
library(mbkmeans)
library(scran)
library(BiocSingular)
library(annotables)
library(tidyverse)
library(SC3)
library(msigdbr)
library(clusterProfiler)
library(flowCore)
library(FlowSOM)
library(factoextra)
library(NbClust)


sce <- read10xCounts("/Volumes/group05/CCBB/CS024892_Kelly_Beshiri/02_PrimaryAnalysisOutput/00_FullCellrangerOutputs/SCAF730_190328_G7/outs/filtered_feature_bc_matrix/")

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
DataFrame(LibSize=sum(qc.lib),
          NExprs=sum(qc.nexprs),
          MitoProp=sum(qc.mito),
          Total=sum(discard))

# Retain only high-quality cells in the SingleCellExperiment.
sce <- sce[,!discard]

## normalization
sce <- normalize(sce)

## feature selection
fit <- trendVar(sce, use.spikes = FALSE)
dec <- decomposeVar(sce, fit)
hvg <- rownames(dec[dec$bio > 0, ]) # ~4k genes

## dimensionality reduction
set.seed(1234)
sce <- runPCA(sce, 
              feature_set = hvg)

sce <- runUMAP(sce)
## SC3 requires this column to be appended
rowData(sce)$feature_symbol <- rowData(sce)$Symbol
## SC3 cannot handle sparse matrixes/hdf nonsense
counts(sce) <- as.matrix(counts(sce))
logcounts(sce) <- as.matrix(logcounts(sce))
# sce <- sc3(sce,
#            ks = 3:6,
#            k_estimator = TRUE)
# 
# sc3_cols <- paste0('sc3_', 3:6, '_clusters')
# 
# p_l <- map(sc3_cols, ~ plotUMAP(sce, colour_by = .))
# 
# wrap_plots(p_l, ncol = 2)
# 
# col_data <- colData(sce)
# head(col_data[ , grep("sc3_", colnames(col_data))])
# 
# plotUMAP(
#   sce, 
#   colour_by = "sc3_5_clusters" 
# )
# sc3_plot_markers(sce, k = 4)
# 
# counts_kmeans <- mbkmeans(sce,
#                           reduceMethod = NA,
#                           whichAssay = "logcounts",
#                           clusters = 10)
# 
# colData(sce)$mbkmeans_3_clusters <- factor(counts_kmeans$Clusters)
# 
# plotUMAP(
#   sce, 
#   colour_by = "mbkmeans_3_clusters" 
# )
# 
# spearman_test <- cor(counts(sce), use = "pairwise.complete", method = "spearman")

logcounts(sce_sc3) <- counts(sce_sc3)
sce_sc3 <- sc3_prepare(sce)
str(metadata(sce_sc3)$sc3)
sce_sc3 <- sc3_estimate_k(sce_sc3)
str(metadata(sce_sc3)$sc3)
sce_sc3 <- sc3_calc_dists(sce_sc3)
names(metadata(sce_sc3)$sc3$distances)

dists_counts <- metadata(sce_sc3)$sc3$distances

dataset <- counts(sce_sc3)
i <- NULL
distances <- c("euclidean", "pearson", "spearman")
message("Calculating distances between the cells...")
n_cores <- 7
cl <- parallel::makeCluster(n_cores, outfile = "")
doParallel::registerDoParallel(cl, cores = n_cores)
dists <- foreach::foreach(i = distances) %dorng% {
  try({
    SC3:::calculate_distance(dataset, i)
  })
}
parallel::stopCluster(cl)
names(dists) <- distances

dist_spearmnan <- as.dist(dists$spearman) 
hclust_spearman <- hclust(dist_spearmnan)

cor_test_res <- lapply(c(seq_along(2:nrow(colData(sce)))),
                       function(col_index){
                         first_col <- col_index - 1
                         cor_vec <- counts(sce)[,first_col]
                         cor_mat <- counts(sce)[,col_index:ncol(counts(sce))]
                         cor_results <- apply(cor_mat, 2, function(x){
                           res <- cor.test(cor_vec,
                                    x,
                                    method = "spearman") %>% 
                             broom::tidy()  
                           return(res)
                         }) %>% 
                           bind_rows()
                        return(cor_results) 
                       })

cor.test(counts(sce)[,1],
         counts(sce)[,2],
         method = "spearman") %>% broom::tidy()

source("single-cell/scrna2019/algs/glmpca.R")

filtered_counts <- counts(sce)
filtered_counts <- filtered_counts[rowSums(filtered_counts) > 0,]

test_glmpca_poi_30 <- glmpca(filtered_counts,
                          30,
                          fam = "poi")

vectorized_loadings <- lapply(test_glmpca_poi_30$loadings, function(dim_x){
  vector_x <- as.vector(dim_x)
  names(vector_x) <- row.names(test_glmpca_poi_30$loadings)
  vector_x <- vector_x[order(vector_x, decreasing = T)]
  return(vector_x)
})

m_df <-msigdbr()
m_t2g = m_df %>% dplyr::select(gs_name, entrez_gene) %>% as.data.frame()

pca_loadings_entrez <- lapply(vectorized_loadings, function(dim_x){
  names(dim_x) <- (data.frame(ensgene = names(dim_x)) %>% 
                     left_join(grch38 %>% 
                                 distinct(ensgene, entrez) %>% 
                                 dplyr::filter(entrez %in% m_t2g$entrez_gene)) %>% 
                     select(entrez) %>% 
                     dplyr::filter(!is.na(entrez)))[[1]]
  dim_x <- dim_x[unique(names(dim_x))]
  dim_x <- sort(dim_x, decreasing = T)
  return(dim_x)
})

pca_loadings_gsea <- lapply(pca_loadings_entrez, GSEA, TERM2GENE = m_t2g, pvalueCutoff = 1)

loadings_gsea_symbols <- lapply(pca_loadings_gsea, function(x){
  new_result <- x@result %>% 
    mutate(entrez = str_split(core_enrichment, "/")) %>%
    unnest(entrez) %>%
    mutate(entrez = as.integer(entrez)) %>% 
    left_join(grch38 %>% 
                dplyr::select(entrez, symbol)) %>% 
    dplyr::select(-c(core_enrichment, entrez)) %>% 
    group_by_at(vars(-symbol)) %>% 
    summarize(core_enrichment = paste(symbol, collapse = ", "))
  return(new_result)
})

gsea_df <- loadings_gsea_symbols %>% 
  bind_rows(.id = "dim")

flow_frame <- new("flowFrame",
                  exprs = as.matrix(test_glmpca_poi_30$factors))

som <- ReadInput(flow_frame)
som <- BuildSOM(som, xdim = 30, ydim = 30)

som_clust <- NbClust(som$map$medianValues, distance = "euclidean", min.nc = 2, max.nc = 10, method = "kmeans")

cluster_mapping <- data.frame(consensus_cluster = som_clust$Best.partition,
                              som_node_id = 1:length(som_clust$Best.partition))

PlotStars(mst)
PlotStars(mst, view = "grid")

save.image("/Volumes/group05/CCBB/CS024892_Kelly_Beshiri/Untitled.RData")

plot_me <- cbind(test_glmpca_poi_30$factors[,1:2], mst$map$mapping[,1]) %>% 
  as.data.frame() %>% 
  setNames(c("dim1", "dim2", "som_node_id")) %>%
  left_join(cluster_mapping)

  ggplot(aes(dim1, dim2, color = factor(consensus_cluster))) +
  geom_point() 

