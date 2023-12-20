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

###############################################################################
#### read in motif enrichments ################################################
###############################################################################

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

###############################################################################
#### read in counts ###########################################################
###############################################################################

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

###############################################################################
#### hnf1a diff bind ##########################################################
###############################################################################

hnf1a_dge_list <- DGEList(counts = hnf1a_counts[,
    str_detect(colnames(hnf1a_counts), ".bam")],
    samples = data.frame(file_name = colnames(
        hnf1a_counts)[str_detect(colnames(hnf1a_counts), ".bam")]) %>%
        mutate(sample_id = str_remove(file_name, "_HNF1A.mLb.clN.sorted.bam"),
        model_id = if_else(str_detect(sample_id, "170"),
          "170",
          "23.1"),
        hnf1a_status = if_else(str_detect(sample_id, "170.2|M"), "lo", "hi"),
            replicate_id = case_when(sample_id %in% c("170.2",
                "H1",
                "M1",
                "170.3T") ~ "r1",
              sample_id %in% c("170.2r2",
                "M1r2",
                "H1r2",
                "170.3r2") ~ "r2",
              T ~ "r3"),
            group_id = paste0(model_id, "_", hnf1a_status)),
            genes = hnf1a_counts[,
                !str_detect(colnames(hnf1a_counts), ".bam")] %>%
                select(Geneid, Length) %>%
                left_join(hnf1a_features,
                    by = c("Geneid" =
                      "PeakID (cmd=annotatePeaks.pl HNF1A.consensus_peaks.bed genome.fa -gid -gtf genes.gtf -cpu 6)")
                    )
    )

hnf1a_dge_list <- hnf1a_dge_list[,
  hnf1a_dge_list$samples$sample_id %in%
    c("M1r2",
      "H2r2",
      "M2r2",
      "H1r2")]

# hnf1a dbgap submission

hnf1a_dge_list$counts %>%
  as.data.frame() %>%
  bind_cols(hnf1a_dge_list$genes[, "Geneid"]) %>%
  setNames(colnames(.) %>%
    str_remove(".mLb.clN.sorted.bam")) %>%
  remove_rownames() %>%
  column_to_rownames("...5") %>%
  write.csv("/data/capaldobj/incoming-nih-dme/CS035088-chip-seq-ilya-hnf1a/hnf1a-counts.csv")

hnf1a_dge_list$samples %>%
    rownames_to_column("RNASEQ.ID") %>%
  transmute(file_id = str_remove(RNASEQ.ID, ".mLb.clN.sorted.bam"),
            model_id = if_else(group == "23.1_hi",
                "HNF1A+ LuCaP 23.1",
                "LuCaP 23.1")) %>%
  write_csv("/data/capaldobj/incoming-nih-dme/CS035088-chip-seq-ilya-hnf1a/hnf1a-file-key.csv")


pdxos <- hnf1a_dge_list[filterByExpr(hnf1a_dge_list,
  group = hnf1a_dge_list$samples$group_id), ]

pdxos <- calcNormFactors(pdxos)

prcomp(t(cpm(pdxos, log = T)))$x %>%
  as.data.frame() %>%
  rownames_to_column("file_name") %>%
  left_join(pdxos$samples) %>%
  ggplot(aes(PC1, PC2)) +
  geom_point() +
  geom_label(aes(label = sample_id, fill = replicate_id))

design <- model.matrix(~0 + hnf1a_status + replicate_id,
                       data = pdxos$samples)
dge_list <- estimateDisp(pdxos, design = design)
fit <- glmQLFit(dge_list, design)

hi_vs_lo_paired <- topTags(glmQLFTest(fit,
  contrast = c(1, -1, 0)),
  n = Inf)$table

###############################################################################
#### h3k27ac diff bind ########################################################
###############################################################################

h3k27ac_dge_list <- DGEList(counts = h3k27ac_counts[,
  str_detect(colnames(h3k27ac_counts), ".bam")],
  samples = data.frame(file_name = colnames(
        h3k27ac_counts)[str_detect(colnames(h3k27ac_counts), ".bam")]) %>%
        mutate(sample_id = str_remove(file_name, "_H3K27Ac.mLb.clN.sorted.bam"),
        model_id = if_else(str_detect(sample_id, "170"),
          "170",
          "23.1"),
        hnf1a_status = if_else(str_detect(sample_id, "170.2|M"), "lo", "hi"),
            replicate_id = case_when(sample_id %in% c("170.2",
                "H1",
                "M1",
                "170.3T") ~ "r1",
              sample_id %in% c("170.2r2",
                "M1r2",
                "H1r2",
                "170.3r2") ~ "r2",
              T ~ "r3"),
            group_id = paste0(model_id, "_", hnf1a_status)),
  genes = h3k27ac_counts[,
                !str_detect(colnames(h3k27ac_counts), ".bam")] %>%
                select(Geneid, Length) %>%
                left_join(h3k27ac_features,
                    by = c("Geneid" =
                      "PeakID (cmd=annotatePeaks.pl H3K27Ac.consensus_peaks.bed genome.fa -gid -gtf genes.gtf -cpu 6)")
                    ))

h3k27ac_dge_list <- h3k27ac_dge_list[,
  !h3k27ac_dge_list$samples$sample_id %in%
    c("H2",
      "M2",
      "H1",
      "M1")]

# hnf1a dbgap submission

h3k27ac_dge_list$counts %>%
  as.data.frame() %>%
  bind_cols(h3k27ac_dge_list$genes[, "Geneid"]) %>%
  setNames(colnames(.) %>%
    str_remove(".mLb.clN.sorted.bam")) %>%
  remove_rownames() %>%
  column_to_rownames("...9") %>%
  write.csv("/data/capaldobj/incoming-nih-dme/CS035088-chip-seq-ilya-hnf1a/h3k27ac-counts.csv")

h3k27ac_dge_list$samples %>%
    rownames_to_column("RNASEQ.ID") %>%
  transmute(file_id = str_remove(RNASEQ.ID, ".mLb.clN.sorted.bam"),
            model_id = case_when(model_id == "23.1" & hnf1a_status == "hi" ~ "HNF1A+ LuCaP 23.1",
              model_id == "23.1" & hnf1a_status != "hi" ~ "LuCaP 23.1",
              model_id != "23.1" & hnf1a_status == "hi" ~ "LuCaP 170.3",
              model_id != "23.1" & hnf1a_status != "hi" ~ "LuCaP 170.2")) %>%
  write_csv("/data/capaldobj/incoming-nih-dme/CS035088-chip-seq-ilya-hnf1a/h3k27ac-file-key.csv")

# fix md5 file

(read_csv("design-hnf1a-chip.csv") %>%
  filter(!str_detect(sample, "^M1_|^H1_|^M2_|^H2_")) %>%
  pivot_longer(starts_with("fastq"), names_to = "fastq_id", values_to = "file_name") %>%
  mutate(file_name = str_remove(file_name,
    "/data/capaldobj/incoming-nih-dme/CS035088-chip-seq-ilya-hnf1a//")))$file_name

read_tsv("md5sum-chip-seq.txt",
  col_names = F) %>%
  separate(X1,
    c("hash",
      "file_name"),
    sep = 34,
    remove = F) %>%
    filter(file_name %in%
      (read_csv("design-hnf1a-chip.csv") %>%
  filter(!str_detect(sample, "^M1_|^H1_|^M2_|^H2_")) %>%
  pivot_longer(starts_with("fastq"), names_to = "fastq_id", values_to = "file_name") %>%
  mutate(file_name = str_remove(file_name,
    "/data/capaldobj/incoming-nih-dme/CS035088-chip-seq-ilya-hnf1a//")))$file_name
) %>% select(X1) %>%
write_tsv("md5sum-chip-seq.txt", col_names = F)

h3k27ac_dge_list <- h3k27ac_dge_list[filterByExpr(h3k27ac_dge_list,
  group = h3k27ac_dge_list$samples$group_id), ]

h3k27ac_dge_list <- calcNormFactors(h3k27ac_dge_list)

prcomp(t(cpm(h3k27ac_dge_list, log = T)))$x %>%
  as.data.frame() %>%
  rownames_to_column("file_name") %>%
  left_join(h3k27ac_dge_list$samples) %>%
  ggplot(aes(PC1, PC2)) +
  geom_point() +
  geom_label(aes(label = sample_id, fill = replicate_id))

h3k27ac_design <- model.matrix(~0 + hnf1a_status + model_id + replicate_id,
                       data = h3k27ac_dge_list$samples)
h3k27ac_dge_list <- estimateDisp(h3k27ac_dge_list, design = h3k27ac_design)
h3k27ac_fit <- glmQLFit(h3k27ac_dge_list, h3k27ac_design)

hi_vs_lo_h3k27ac_paired <- topTags(glmQLFTest(h3k27ac_fit,
  contrast = c(1, -1, 0, 0, 0)),
  n = Inf)$table

h3k27ac_design_split <- model.matrix(~0 + group_id + replicate_id,
                       data = h3k27ac_dge_list$samples)
h3k27ac_dge_list_split <- estimateDisp(h3k27ac_dge_list,
  design = h3k27ac_design_split)
h3k27ac_fit_split <- glmQLFit(h3k27ac_dge_list_split,
  h3k27ac_design_split)

hi_vs_lo_h3k27ac_170 <- topTags(glmQLFTest(h3k27ac_fit_split,
  contrast = c(1, -1, 0, 0, 0, 0)),
  n = Inf)$table

hi_vs_lo_h3k27ac_23.1 <- topTags(glmQLFTest(h3k27ac_fit_split,
  contrast = c(0, 0, 1, -1, 0, 0)),
  n = Inf)$table

###############################################################################
#### plots ####################################################################
###############################################################################

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