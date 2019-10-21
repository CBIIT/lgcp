library(monocle)
library(scran)
library(scater)

sce_adeno <- sce_glm_pca[,reducedDim(sce_glm_pca, "UMAP")[,1] > -5]

cds <- convertTo(sce_adeno, type = "monocle",
                 row.fields = c(1:ncol(rowData(sce_adeno))),
                 col.fields = c(1:ncol(colData(sce_adeno))),
                 expressionFamily = negbinomial.size())

featureData(cds)$gene_short_name <- rowData(sce_adeno)$Symbol

cds <- estimateSizeFactors(cds) 
cds <- estimateDispersions(cds) 
cds <- detectGenes(cds)

cds <- setOrderingFilter(cds, row.names(fData(cds))[fData(cds)$dev <= 2000])
cds <- reduceDimension(cds, method = "DDRTree")
cds <- orderCells(cds)

colData(sce_adeno) <- cbind(colData(sce_adeno), pData(cds)$Pseudotime, pData(cds)$State)

plotReducedDim(sce_adeno, "UMAP", colour_by = "pData(cds)$Pseudotime")
