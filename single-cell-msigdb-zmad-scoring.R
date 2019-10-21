SCAF730_190328_G7 = MB44 (G7)
SCAF731_190328_G18 = LuCaP 145.2 (G18)
SCAF732_190328_G6 = LuCaP 173.1 (G6)

library(msigdbr)

logcounts(sce_glm_pca) %>% 
  tidy_sparse_matrix() %>% 
  left_join(counts(sce_glm_pca) %>% 
              tidy_sparse_matrix())

# calculate mean absolute deviation for each gene across all samples
median_absolute_deviations <- logcounts(sce_glm_pca) %>% 
  tidy_sparse_matrix() %>%
  group_by(row) %>%
  summarise(median_abs_dev = mad(value),
            median = median(value))

# calculate modified Z-score by multipling the difference of the cpm_tmm and mean for each gene times 0.6745 (just because), and dividing by the mean absolute deviation
modified_z_score <- logcounts(sce_glm_pca) %>% 
  tidy_sparse_matrix() %>%
  left_join(median_absolute_deviations) %>%
  mutate(mod_z = (0.6745*(value - median))/median_abs_dev) %>% 
  select(row, column, mod_z)

fake_cell <- as.integer(ncol(sce_glm_pca) + 1)

zmad_sparse <- data.frame(row = names(logcounts(sce_glm_pca)[,1])) %>% 
  left_join(modified_z_score) %>% 
  mutate(row = factor(row, levels = names(logcounts(sce_glm_pca)[,1])),
         column = if_else(is.na(column), fake_cell, column),
         mod_z = if_else(is.na(mod_z), 0, mod_z)) %>%
  arrange(column) %>% 
  tidytext::cast_sparse(row = row,
                        column = column,
                        value = mod_z)

dimnames(zmad_sparse) <- list(dimnames(zmad_sparse)[[1]], NULL)
assay(sce_glm_pca, "zmad") <- zmad_sparse[,-fake_cell]

# create function to calculate scores
sc_scores <- msigdbr() %>%
  left_join(modified_z_score,
            by = c("entrez_gene" = "entrez")) %>%
  group_by(gs_name, column) %>%
  summarise(ont_score = mean(mod_z))

sc_score <- lapply(seq_along(1:ncol(sce_glm_pca)), function(i){
  scored <- data.frame(mod_z = assay(sce_glm_pca, "zmad")[,i],
                       ensgene = names(assay(sce_glm_pca, "zmad")[,i])) %>%
    filter(mod_z > 0) %>% 
    left_join(grch38 %>% 
                select(ensgene, entrez) %>% 
                distinct()) %>%
    left_join(msigdbr(),
              by = c("entrez" = "entrez_gene")) %>% 
    group_by(gs_name) %>% 
    summarize(score = mean(mod_z))
  return(scored)
})
