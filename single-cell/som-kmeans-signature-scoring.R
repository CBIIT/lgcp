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
                             })

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
                      })

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
                         })

barcode_ar_signature_score <- vector(length = ncol(sce_glm_pca))
for(i in 1:ncol(sce_glm_pca)){
  scored <- data.frame(logcounts = logcounts(sce_glm_pca)[rownames(logcounts(sce_glm_pca)) %in% ar_signature_weights$`Ensemble ID`,i],
                       ensgene = names(logcounts(sce_glm_pca)[rownames(logcounts(sce_glm_pca)) %in% ar_signature_weights$`Ensemble ID`,i])) %>% 
    left_join(ar_signature_weights, by = c("ensgene" = "Ensemble ID")) %>% 
    transmute(score = logcounts*`BinReg Coef`)
  
  barcode_ar_signature_score[i] <- sum(scored$score)
}
# NE score
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

# AR score
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
