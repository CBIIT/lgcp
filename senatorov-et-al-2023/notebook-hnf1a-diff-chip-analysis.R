#' ---
#' title: "Bone mets project single cell RNA sequencing compilation data analysis"
#' author: "Brian Capaldo"
#' date: "10/06/2022"
#' output:
#'   html_document:
#'     toc: true
#'     toc_float: false
#'     number_sections: true
#'     code_folding: hide
#' ---

#+ message = FALSE, warning = FALSE

###############################################################################
#### load up libraries ########################################################
###############################################################################

library(tidyverse)
library(edgeR)
library(annotables)
library(readxl)
library(DT)
library(SummarizedExperiment)
library(gplots)
library(pvclust)
library(clusterProfiler)
library(msigdbr)
library(GSVA)
library(variancePartition)
library(limma)
library(BiocParallel)
library(gplots)

param <- SnowParam(8, "SOCK", progressbar = FALSE)
register(param)

set.seed(8675309)

motif_matches <- list.files(
    paste0("/data/capaldobj/incoming-nih-dme/CS035088-chip-seq-ilya-hnf1a/",
        "chipseq-grch37-results/MotifMatching/"),
           full.names = T) %>%
  lapply(read_tsv, col_names = c("chr",
    "start",
    "end",
    "motif_id",
    "score",
    "strand"))

motif_matches_list <- motif_matches %>%
  setNames(list.files(
    paste0("/data/capaldobj/incoming-nih-dme/CS035088-chip-seq-ilya-hnf1a/",
        "chipseq-grch37-results/MotifMatching/")
  )) %>%
  bind_rows(.id = "sample_id") %>%
  pivot_wider(names_from = sample_id,
              values_from = score) %>%
  split(.$motif_id)

hnf1a_counts <- read_tsv(
    paste0("/data/capaldobj/incoming-nih-dme/CS035088-chip-seq-ilya-hnf1a/",
        "chipseq-grch37-results/bwa/mergedLibrary/macs2/narrowPeak/consensus/",
        "HNF1A/HNF1A.consensus_peaks.featureCounts.txt"),
         skip = 1)

hnf1a_features <- read_tsv(
    paste0("/data/capaldobj/incoming-nih-dme/CS035088-chip-seq-ilya-hnf1a/",
        "chipseq-grch37-results/bwa/mergedLibrary/macs2/narrowPeak/consensus/",
        "HNF1A/HNF1A.consensus_peaks.annotatePeaks.txt"))

h3k27ac_counts <- read_tsv(
    paste0("/data/capaldobj/incoming-nih-dme/CS035088-chip-seq-ilya-hnf1a/",
        "chipseq-grch37-results/bwa/mergedLibrary/macs2/narrowPeak/consensus/",
        "H3K27Ac/H3K27Ac.consensus_peaks.featureCounts.txt"),
                           skip = 1)

h3k27ac_features <- read_tsv(
    paste0("/data/capaldobj/incoming-nih-dme/CS035088-chip-seq-ilya-hnf1a/",
        "chipseq-grch37-results/bwa/mergedLibrary/macs2/narrowPeak/consensus/",
        "H3K27Ac/H3K27Ac.consensus_peaks.annotatePeaks.txt"))

hnf1a_dge_list <- DGEList(counts = hnf1a_counts[,
    str_detect(colnames(hnf1a_counts), ".bam")],
    samples = data.frame(file_name = colnames(
        hnf1a_counts)[str_detect(colnames(hnf1a_counts), ".bam")]) %>%
        mutate(sample_id = str_remove(file_name, "_HNF1A.mLb.clN.sorted.bam"),
        model_id = c("170",
            "23.1",
            "23.1",
            "170",
            "23.1",
            "23.1",
            "170",
            "170",
            "23.1",
            "23.1"),
        hnf1a_status = c("lo",
            "hi",
            "lo",
            "lo",
            "hi",
            "hi",
            "hi",
            "hi",
            "lo",
            "lo"),
            replicate_id = paste0("r", c(1,1,2,2,3,2,2,1,3,1)), 
            group_id = paste0(model_id, "_", hnf1a_status)),
            genes = hnf1a_counts[,
                !str_detect(colnames(hnf1a_counts), ".bam")] %>%
                select(Geneid, Length) %>%
                left_join(hnf1a_features,
                    by = c("Geneid" = "PeakID (cmd=annotatePeaks.pl HNF1A.consensus_peaks.bed genome.fa -gid -gtf genes.gtf -cpu 6)")
                    )
    )

pdxos <- hnf1a_dge_list[filterByExpr(hnf1a_dge_list, group = hnf1a_dge_list$samples$group_id),]
pdxos <- calcNormFactors(pdxos)

prcomp(t(cpm(pdxos, log = T)))$x %>% 
  as.data.frame() %>% 
  rownames_to_column("file_name") %>% 
  left_join(pdxos$samples) %>% 
  ggplot(aes(PC1, PC2)) +
  geom_point() +
  geom_label(aes(label = sample_id, fill = replicate_id))

design <- model.matrix(~0 + hnf1a_status + model_id + replicate_id,
                       data = pdxos$samples)
dge_list <- estimateDisp(pdxos, design = design)
fit <- glmQLFit(dge_list, design)

hi_vs_lo_paired <- topTags(glmQLFTest(fit, contrast = c(1,-1,0,0,0)), n = Inf)$table 

hi_vs_lo_paired %>% 
  ggplot(aes(FDR)) +
  geom_histogram(binwidth = 0.05)

ggplot(hi_vs_lo_paired, aes(logFC, -log10(FDR))) +
  geom_point()

cpm(pdxos, log = T) %>%
  as.data.frame() %>%
  bind_cols(pdxos$genes) %>%
  filter(Geneid == "Interval_13313") %>%
  pivot_longer(ends_with("bam"),
               names_to = "file_name",
               values_to = "expression") %>%
  left_join(pdxos$samples) %>%
  ggplot(aes(group_id, expression)) +
  geom_point() +
  ggrepel::geom_label_repel(aes(label = replicate_id))

cpm(pdxos, log = T) %>%
  as.data.frame() %>%
  bind_cols(pdxos$genes) %>%
  filter(Geneid == "Interval_54591") %>%
  pivot_longer(ends_with("bam"),
               names_to = "file_name",
               values_to = "expression") %>%
  left_join(pdxos$samples) %>%
  ggplot(aes(group_id, expression)) +
  geom_point() +
  ggrepel::geom_label_repel(aes(label = replicate_id))

design <- model.matrix(~0 + group_id + replicate_id,
                       data = pdxos$samples)
dge_list <- estimateDisp(pdxos, design = design)
fit <- glmQLFit(dge_list, design)

hi_vs_lo_170 <- topTags(glmQLFTest(fit, contrast = c(1,-1,0,0,0,0)), n = Inf)$table 
hi_vs_lo_23.1 <- topTags(glmQLFTest(fit, contrast = c(0,0,1,-1,0,0)), n = Inf)$table 

list("paired" = hi_vs_lo_paired,
     "lucap_23.1" = hi_vs_lo_23.1,
     "lucap_170" = hi_vs_lo_170) %>% 
  bind_rows(.id = "model_id") %>% 
  ggplot(aes(FDR)) +
  geom_histogram(binwidth = 0.05) +
  facet_wrap(~model_id)

list("paired" = hi_vs_lo_paired,
     "lucap_23.1" = hi_vs_lo_23.1,
     "lucap_170" = hi_vs_lo_170) %>% 
  bind_rows(.id = "model_id") %>% 
  ggplot(aes(logFC, -log10(FDR))) +
  geom_point() +
  facet_wrap(~model_id)
