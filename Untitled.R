## feature selection
fit <- trendVar(sce, use.spikes = F)
dec <- decomposeVar(sce, fit)
hvg <- rownames(dec[dec$bio > 0, ]) # ~4k genes

## dimensionality reduction
sce_pca <- runPCA(sce, 
                  feature_set = hvg,
                  ncomponents = 30)

sce_pca <- runUMAP(sce_pca)

g_pca <- buildSNNGraph(sce_pca, k=10, use.dimred = 'PCA')
clust_pca <- igraph::cluster_walktrap(g_pca)$membership
table(clust_pca)
sce_pca$clust <- factor(clust_pca)
sce_pca$NE_score <- barcode_neuro_score
sce_pca$AR_score <- barcode_ar_signature_score


plotReducedDim(sce_glm_pca, "UMAP", colour_by = "AR_score")
plotReducedDim(sce_glm_pca, "UMAP", colour_by = "NE_score")
plotReducedDim(sce_glm_pca, "UMAP", colour_by = "clust")

plotReducedDim(sce_pca, "PCA", ncomponents = 5, colour_by = "clust")
plotReducedDim(sce_glm_pca, "GLM_PCA", ncomponents = 5, colour_by = "clust")

plotReducedDim(sce_glm_pca, "GLM_PCA")


load("~/Desktop/single-cell/ct-35-2-v1-sce-object.Rdata")
sce_glm_pca_ct35_2 <- sce_glm_pca

load("~/Desktop/single-cell/lucap173.1-sce-object.Rdata")
sce_glm_pca_lucap_173.1 <- sce_glm_pca

load("~/Desktop/single-cell/mb155-sce-object.Rdata")

plotReducedDim(sce_glm_pca_ct35_2, "GLM_PCA", colour_by = c("som_dim1_dim2")) + theme(legend.position = "none")
plotReducedDim(sce_glm_pca_ct35_2, "GLM_PCA", colour_by = c("som_dim1_dim2")) + theme(legend.position = "none")

plotReducedDim(sce_glm_pca_ct35_2, "UMAP", colour_by = c("som_dim1_dim2")) + theme(legend.position = "none")

plotReducedDim(sce_glm_pca_ct35_2, "GLM_PCA", colour_by = c("AR_score"), ncomponents = 2:3)
plotReducedDim(sce_glm_pca_ct35_2, "GLM_PCA", colour_by = c("NE_score"), ncomponents = 1:2)

CT35_2 <- colData(sce_glm_pca_ct35_2) %>% 
  as.data.frame() %>% 
  select(AR_score, NE_score)

LuCaP_173.1 <- colData(sce_glm_pca_lucap_173.1) %>% 
  as.data.frame() %>% 
  select(AR_score, NE_score)

MB155 <- colData(sce_glm_pca) %>% 
  as.data.frame() %>% 
  select(AR_score, NE_score)

list("CT35-2" = CT35_2,
     "LuCaP 173.1" = LuCaP_173.1,
     "MB155" = MB155) %>% 
  bind_rows(.id = "sample_id") %>% 
  ggplot(aes(AR_score, NE_score)) +
  geom_point() +
  facet_grid(~sample_id) +
  theme_bw()

plotReducedDim(sce_glm_pca, "UMAP", colour_by = "clust")
plotReducedDim(sce_glm_pca_lucap_173.1, "UMAP", colour_by = "clust")
f  
