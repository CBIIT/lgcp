library(DropletUtils)
library(scater)
library(scran)
library(BiocSingular)
library(annotables)
library(tidyverse)
library(msigdbr)
library(clusterProfiler)
library(flowCore)
library(FlowSOM)
source("single-cell/scrna2019/algs/glmpca.R")


# read in 10x counts, change this to parameter fed in at command line
#sce <- read10xCounts("/Volumes/group05/CCBB/CS024892_Kelly_Beshiri/02_PrimaryAnalysisOutput/00_FullCellrangerOutputs/SCAF730_190328_G7/outs/filtered_feature_bc_matrix/")

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
## SC3 requires this column to be appended
rowData(sce)$feature_symbol <- rowData(sce)$Symbol
## SC3 cannot handle sparse matrixes/hdf nonsense
counts(sce) <- as.matrix(counts(sce))
logcounts(sce) <- as.matrix(logcounts(sce))

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
    dplyr::distinct() %>% 
    group_by_at(vars(-symbol)) %>% 
    summarize(num_genes = n_distinct(symbol),
      core_enrichment = paste(symbol, collapse = ", "))
  return(new_result)
})


gsea_df <- loadings_gsea_symbols %>% 
  bind_rows(.id = "dim")


test_clust <- lapply(test_glmpca_poi_30$factors, kmeans, centers = 2, iter.max = 100, nstart = 10)

plot_df <- data.frame(test_glmpca_poi_30$factors,
                      Barcode = colData(sce)$Barcode) %>% 
  gather(reduced_dim_id, value, -Barcode) %>% 
  bind_cols(lapply(test_clust, `[[`, 1) %>%
              lapply(as.data.frame) %>% 
              bind_rows(.id = "dim_id") %>% 
              setNames(c("dim_id", "kmeans_id"))) %>% 
  mutate(reduced_dim_id = factor(reduced_dim_id, levels = paste0("dim", c(1:30))))

kmeans_dim_id <- plot_df %>%
  group_by(reduced_dim_id, kmeans_id) %>%
  summarize(cluster_median = median(value)) %>%
  arrange(reduced_dim_id, cluster_median) %>%
  ungroup() %>% 
  mutate(cluster_id = factor(rep(c("low", "high"), 30), levels = c("low", "high")))

plot_df <- plot_df %>% 
  left_join(kmeans_dim_id) 

plot_df %>% 
  ggplot(aes(value, fill = cluster_id)) +
  geom_histogram(bins = 100) +
  facet_wrap(~reduced_dim_id) +
  theme_bw()

flow_frame <- new("flowFrame",
                  exprs = as.matrix(test_glmpca_poi_30$factors))

som <- ReadInput(flow_frame)
som <- BuildSOM(som, xdim = 30, ydim = 30)

plot_df <- plot_df %>% 
  left_join(data.frame(Barcode = colData(sce)$Barcode,
                       som_node = mst$map$mapping[,1])) %>% 
  left_join(data.frame(som_node = c(1:nrow(mst$MST$l)),
            mst_x = mst$MST$l[,1],
            mst_y = mst$MST$l[,2]))

plot_df %>%
  ggplot(aes(mst_x, mst_y, color = cluster_id)) +
  geom_point() +
  facet_wrap(~reduced_dim_id) +
  theme_void() +
  theme(legend.position = "none")

som_nodes_stats <- plot_df %>% 
  group_by(som_node, reduced_dim_id) %>% 
  summarise(som_mean = mean(value),
            som_median = median(value)) %>%
  left_join(plot_df %>% 
              group_by(som_node, reduced_dim_id) %>%
              count(cluster_id) %>% 
              ungroup() %>% 
              spread(cluster_id,
                     n,
                     fill = 0) %>% 
              mutate(fraction_high = high/(low+high),
                     num_cells = low+high)) %>% 
  left_join(kmeans_dim_id %>% 
              dplyr::select(reduced_dim_id, cluster_id, cluster_median) %>% 
              spread(cluster_id, cluster_median) %>% 
              setNames(c("reduced_dim_id", "kmns_lo_clstr_med", "kmns_hi_clstr_med")))

playing <- som_nodes_stats %>% 
  mutate(som_clst_label = case_when(som_median <= kmns_lo_clstr_med ~ "Negative",
                                    som_median >= kmns_hi_clstr_med ~ "Positive",
                                    som_median > kmns_lo_clstr_med && fraction_high < 0.25 ~ "Low",
                                    som_median < kmns_hi_clstr_med && fraction_high > 0.75 ~ "High",
                                    T ~ "Neutral" ),
         som_clst_label = factor(som_clst_label, levels = c("Negative", "Low", "Neutral", "High", "Positive")))

playing %>%
  left_join(data.frame(som_node = c(1:nrow(mst$MST$l)),
                       mst_x = mst$MST$l[,1],
                       mst_y = mst$MST$l[,2])) %>% 
  ggplot(aes(mst_x, mst_y, color = som_clst_label, size = num_cells, alpha = 0.25)) +
  geom_point(shape = 21) +
  facet_wrap(~reduced_dim_id) +
  theme_void() +
  scale_color_manual(values = c("darkblue", "lightblue", "grey", "red", "darkred")) +
  theme(legend.position = "none")

ggsave("/Volumes/group05/CCBB/CS024892_Kelly_Beshiri/som-plot.pdf",
       width = 16,
       height = 9)

gsea_df %>% 
  write_csv("/Volumes/group05/CCBB/CS024892_Kelly_Beshiri/glm-pcs-gsea-results.csv")

playing %>% 
  dplyr::select(reduced_dim_id, som_clst_label, som_node, num_cells) %>% 
  spread(reduced_dim_id, som_clst_label) %>% 
  ungroup() %>%
  arrange(desc(num_cells))
  dplyr::select(-som_node) %>% 
  distinct()

save.image("/Volumes/group05/CCBB/CS024892_Kelly_Beshiri/Untitled.RData")

plot_me <- cbind(test_glmpca_poi_30$factors[,1:2], mst$map$mapping[,1]) %>% 
  as.data.frame() %>% 
  setNames(c("dim1", "dim2", "som_node_id")) %>%
  left_join(cluster_mapping)

  ggplot(aes(dim1, dim2, color = factor(consensus_cluster))) +
  geom_point() 

