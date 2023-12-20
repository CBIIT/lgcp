#' ---
#' title: "23.1 HNF analysis"
#' author: "Brian Capaldo"
#' date: "02/25/2022"
#' output:
#'   html_document:
#'     toc: true
#'     toc_float: true
#'     number_sections: true
#'     code_folding: hide
#' ---

#+ message = FALSE, warning = FALSE

library(tidyverse)
library(edgeR)
library(annotables)
library(readxl)
library(msigdbr)
library(clusterProfiler)
library(DT)
library(GSVA)

################################################################################
#### load data #################################################################
################################################################################

#' # Results
#'
#' Quantified counts generated using the nextflow core RNA sequencing pipeline.
#' Briefly, fastq files were aligned against both grch37 using STAR and reads
#' were quantified using featureCounts from Subread. Raw fastq files were
#' trimmed using TrimGalore!. Low expressed genes have been filtered out using
#' `edgeR::filterByExpr`. Basically, a gene is kept if a meaningful number of
#' samples within each group (in this case, groups are based condition) show
#' expression of the gene. A gene must be expressed in nearly all replicates to
#' be kept. This does mean some genes will be kept that have no expression
#' in some groups.

#+ message = F, warning = F

# read in counts
grch37_counts <- read_tsv(
    "/data/capaldobj/incoming-nih-dme//CS031135-drug-screen-ko-lncap-crispr-ko-deg/human-samples-hg19-rna-seq-results//featureCounts/merged_gene_counts.txt"
)
# grab 23.1 samples only
counts <- grch37_counts[,c(1:2, grep("IS0", colnames(grch37_counts)))]
# finish building dge_list
dge_list_grch37 <- DGEList(counts = as.matrix(counts[,-c(1:2)]),
    genes = counts[,c(1:2)] %>%
        left_join(grch37 %>%
        distinct(ensgene, .keep_all = T),
            by = c("Geneid" = "ensgene")))

dge_list_grch37$samples <- dge_list_grch37$samples %>% 
  rownames_to_column("file_name") %>%
  mutate(sample_id = str_remove(file_name, "_S.*"),
    sample_condition = if_else(sample_id %in% c("IS01", "IS02", "IS03", "IS04"),
        "KI", "EV") %>%
        factor(levels = c("EV", "KI"))) %>%
  column_to_rownames("file_name")

# output for dbgap submission

dge_list_grch37$counts %>%
  as.data.frame() %>%
  bind_cols(dge_list_grch37$genes[, "Geneid"]) %>%
  setNames(colnames(.) %>%
    str_remove("Aligned.sortedByCoord.out.bam")) %>%
  remove_rownames() %>%
  column_to_rownames("...9") %>%
  write.csv("/data/capaldobj/incoming-nih-dme//CS031135-drug-screen-ko-lncap-crispr-ko-deg/grch37-rna-counts.csv")

dge_list_grch37$samples %>%
    rownames_to_column("RNASEQ.ID") %>%
  transmute(file_id = str_remove(RNASEQ.ID, "Aligned.sortedByCoord.out.bam"),
            model_id = if_else(sample_condition == "KI",
                "HNF1A+ LuCaP 23.1",
                "LuCaP 23.1")) %>%
  write_csv("/data/capaldobj/incoming-nih-dme//CS031135-drug-screen-ko-lncap-crispr-ko-deg/file-key.csv")

dge_list_grch37 <- calcNormFactors(dge_list_grch37[
    filterByExpr(dge_list_grch37, group = dge_list_grch37$samples$sample_condition),
])

################################################################################
#### PCA #######################################################################
################################################################################

#' ## Principal components analysis
#' 
#' Principal components analysis (PCA) is a common approach for evaluating the
#' degree of similarity between samples using variance of the expression across
#' samples. In this case, we run PCA on the human and mouse reads separately.
#' The end goal is to crete a means of modeling the tumor microenvironment as a
#' series of principal component embeddings which represent the similarity of
#' the transcriptomes between samples.

#+ message = F, warning = F
human_pca <- prcomp(t(cpm(dge_list_grch37, log = T)))

#' ### Human PCA plot

#+ message = F, warning = F

dge_list_grch37$samples %>% 
    rownames_to_column("file_name") %>%
    left_join(human_pca$x %>% 
        as.data.frame() %>% 
        rownames_to_column("file_name") %>% 
        select(file_name, starts_with("PC"))) %>%
  ggplot(aes(PC1, PC2)) +
  geom_point() +
  ggrepel::geom_label_repel(aes(label = sample_id)) +
  theme_bw()

################################################################################
#### DEG #######################################################################
################################################################################

#' ## Differential expression analysis
#'
#' Differential expression analysis was performed to compare the knock-in (KI)
#' against the empty vector (EV). As expected, amongst the most upregulated in
#' the KI was HNF1A. It should be noted that the KI was quite severe in terms of
#' its effects. Looking at the histogram of FDRs, we can see nearly 7500 genes
#' showing DGE with an FDR of better than 0.05 (the first column in the
#' histogram).

#+ message = F, warning = F

design <- model.matrix(~0 + sample_condition,
                       data = dge_list_grch37$samples)

dge_list <- estimateDisp(dge_list_grch37, design = design)
fit <- glmQLFit(dge_list, design)

topTable <- topTags(glmQLFTest(fit, contrast = c(-1,1)), n = Inf)$table

topTable %>% 
  ggplot(aes(FDR)) + 
  geom_histogram(binwidth = 0.05)

#' ### Volcano plot
#' 
#' Below is two volcano plots showing a summary of the DGE analysis. Here we
#' are showing log2 fold change against -log10 of the FDR, and the top 20 most
#' upregulated and downregulated genes have been labeled in the first plot.
#' In the second plot, we've labelled genes of interest from the UGT, GLY, and
#' HNF1A. A red line has been drawn showing genes exceeding a FDR of better than
#' 5%. It should be noted that generally, you observe no expression of UGT
#' family members in the EV, to pronounced expression in the KI.

#+ message = F, warning = F, fig.width = 16, fig.height = 9 

topTable %>% 
  arrange(logFC) %>%
  mutate(row_number = row_number(),
    gene_label = if_else(row_number <= 20 | row_number >= 21829, symbol, "")) %>%
  ggplot(aes(logFC, -log10(FDR))) + 
  geom_point() +
  ggrepel::geom_label_repel(aes(label = gene_label), max.overlaps = 100) +
  geom_hline(yintercept = -log10(0.05), color = "red", alpha = 0.5)

topTable %>% 
  mutate(gene_label = if_else(str_detect(symbol, "UGT|GLY|HNF1A"), symbol, "")) %>%
  ggplot(aes(logFC, -log10(FDR))) + 
  geom_point() +
  ggrepel::geom_label_repel(aes(label = gene_label), max.overlaps = 1000) +
  geom_hline(yintercept = -log10(0.05), color = "red", alpha = 0.5)

################################################################################
#### PEA #######################################################################
################################################################################

msigdb <- msigdbr()

run_gsea_plots <- function(sig_top_table, msigdb){
  sig_top_table <- sig_top_table
  enrichment_contrasts <- lapply(unique(sig_top_table$contrast_id), function(contrast){
    sig_cont_top_table <- sig_top_table[sig_top_table$contrast_id == contrast,]
    enrichment_cats <- lapply(unique(msigdb$gs_cat), function(category){
      gsea_gene_list <- sig_cont_top_table %>%
        arrange(desc(logFC))

      gsea_gene_list <- setNames(gsea_gene_list$logFC, gsea_gene_list$symbol)

      gseMSIGDB_results <- GSEA(gsea_gene_list,
                                TERM2GENE = msigdb %>%
                                  dplyr::filter(gs_cat == category) %>%
                                  transmute(ont = gs_name,
                                            gene = gene_symbol))

      return(gseMSIGDB_results)
    })
    names(enrichment_cats) <- unique(msigdb$gs_cat)
    return(enrichment_cats)
  })
  names(enrichment_contrasts) <- unique(sig_top_table$contrast_id)
  return(enrichment_contrasts)
}

#' ### DEG Table
#' 
#' Below is the table of differentially expressed genes.
#' 
#' The column descriptions are:
#' 
#' | column name | description |
#' |:------------|:------------|
#' | Geneid | ensembl human id |
#' | gene_name | gene name |
#' | entrez | entrez identifier |
#' | chr | chromosome |
#' | start | start location of transcript |
#' | end | end location of transcript |
#' | strand | strandedness |
#' | biotype | transcript type |
#' | description | gene description |
#' | logFC | log 2 fold change for the contrast |
#' | logCPM | average expression of the gene across all samples (log2 CPM) |
#' | F | value of the quasiliklihood F statistic (analogous to a t statistic) |
#' | PValue | p value of the contrast |
#' | FDR | false discovery rate adjusted p value for the contrast |
#' 

#+ message = F, warning = F

topTable %>% 
  select(-symbol) %>%
  datatable(extensions = 'Buttons', options = list(
    dom = 'Bfrtip',
    buttons = 
      list('copy', 'print', list(
        extend = 'collection',
        buttons = c('csv', 'excel', 'pdf'),
        text = 'Download'
      ))
    
  )) %>% 
  formatRound(c('logFC', 'logCPM', 'F')) %>% 
  formatSignif(c('PValue', 'FDR'), digits = 4)

#+ warning = F, message = F, fig.width = 16, fig.height = 9

#' ### Gene set enrichment analysis
#' 
#' With the set of all genes, we can also run a gene set enrichment analysis
#' (GSEA). GSEA takes a rank ordered list of genes, and tests if any gene sets 
#' members are not randomly distributed through the rank ordered list. Gene sets
#' with nonrandom orders will have more than expected genes at the front or back 
#' of the list, and the positioning of the genes will be reflected by the 
#' scores; a positive score indicates the gene set has lots of members with high
#' fold changes, putting it at the front of the list, and a negative score 
#' indicates the opposite. The enrichment was performed on every category in the 
#' molecular signatures database (MSigDB) independently. So each q value and 
#' adjusted p value are corrected for the size of the category, and not the 
#' entirety of MSigDB.
#' 
#' | column name | description |
#' |:------------|:------------|
#' | contrast_id | the contrast being made, use table above for details of samples being compared |
#' | geneset_collections_id | [MSigDB gene set collection](https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp) |
#' | ID | gene set name |
#' | Description | link to gene set description |
#' | setSize | number of genes in the gene set |
#' | enrichmentScore | GSEA enrichment score |
#' | NES | normalized GSEA enrichment score |
#' | pvalue | pvalue for the enrichment test |
#' | p.adjust | Benjamani Hochberg (FDR) corrected p value |
#' | qvalues | q adjusted p value |
#' | rank | rank of the gene set |
#' | leading_edge | leading edge statistics |  
#' | core_enrichment | list of well ranked genes present in the gene set |
#' | gs_subcat | [subcategory of gene set](https://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp) |
#'

#+ message = F, warning = F

gsea_results <- run_gsea_plots(topTable %>%
    mutate(contrast_id = "KI vs EV"), msigdb)

gsea_plots_table <- lapply(gsea_results, function(contrast_gsea_results){
  gsea_results <- lapply(contrast_gsea_results, function(enrich_results){
    return(enrich_results@result)
  }) %>% bind_rows(.id = "geneset_collections_id")
}) %>% bind_rows(.id = "contrast_id") %>% 
  left_join(msigdb %>%
              dplyr::distinct(gs_subcat, gs_name),
            by = c("ID" = "gs_name"))

gsea_plots_table %>%
  mutate(Description = paste0("<a href='https://www.gsea-msigdb.org/gsea/msigdb/cards/",
                              Description,
                              "'>",
                              Description,
                              "</a>")) %>%
  datatable(escape = F,
            extensions = 'Buttons', options = list(
              dom = 'Bfrtip',
              buttons =
                list('copy', 'print', list(
                  extend = 'collection',
                  buttons = c('csv', 'excel', 'pdf'),
                  text = 'Download'
                ))
            )) %>%
  formatRound(c('enrichmentScore', 'NES')) %>%
  formatSignif(c('pvalue', 'p.adjust', 'qvalues'), digits = 4)

#+ message = F, warning = F, fig.width = 16, fig.height = 9 

cnetplot(gsea_results$`KI vs EV`$H,
    foldChange = gsea_results$`KI vs EV`$H@geneList) +
    scale_fill_gradient2(low = "blue", high = "red")

cnetplot(gsea_results$`KI vs EV`$C2,
    foldChange = gsea_results$`KI vs EV`$H@geneList) +
    scale_fill_gradient2(low = "blue", high = "red")

cnetplot(gsea_results$`KI vs EV`$C3,
    foldChange = gsea_results$`KI vs EV`$H@geneList) +
    scale_fill_gradient2(low = "blue", high = "red")

save.image("HNF1A-knock-in-dge-workspace.RData")