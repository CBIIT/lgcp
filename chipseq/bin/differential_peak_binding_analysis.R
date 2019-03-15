#!/usr/bin/env Rscript

require(docopt)
require(methods)

"
Usage:
differential_peak_binding_analysis.R (-h | --help | --version)
differential_peak_binding_analysis.R MACS_SAMPLE_SHEET

Description:   This program is a command line interface to edgeR

Options:

--version                  Show the current version.

Arguments:

MACS_SAMPLE_SHEET   Provide macs sample sheet for dba

" -> doc

args <- docopt(doc)

# load up libraries

library(DiffBind)
library(tidyverse)
library(readxl)
library(ChIPpeakAnno)
library(EnsDb.Hsapiens.v75)

annoData <- toGRanges(EnsDb.Hsapiens.v75, feature="gene")

# list all bam files in the picard sub directory
bam_files <- list.files(path = "combined-chip-seq/picard/",
                        pattern = ".bam$",
                        full.names = T)

# read in sample sheet from Nichelle (awkward formating)
sample_sheet <- read_xlsx("ChIP sample sheet.xlsx") %>%
  slice(3:17) %>%
  select(1:3) %>%
  setNames(c("sample_id", "condition", "antibody"))

# turn bam file vector into data frame
bam_files <- data.frame(bamReads = bam_files) %>%
  mutate(sample_id = basename(as.character(bamReads)) %>%
           str_remove("\\.dedup\\.sorted\\.bam"))

peak_files <- data.frame(Peaks = list.files("combined-chip-seq/macs/",
                                            pattern = ".xls",
                                            full.names = T)) %>%
  mutate(sample_id = basename(as.character(Peaks)) %>%
           gsub("sample\\_", "A2-167-", .) %>%
           str_remove("_peaks.xls"))

# join list of bam files with experimental data
sample_sheet <- sample_sheet %>%
  left_join(bam_files) %>%
  left_join(peak_files %>%
              mutate(PeakCaller = rep("macs", nrow(.))))

sample_sheet %>%
  setNames(c("SampleID", "Condition", "Factor", "bamReads", "Peaks", "PeakCaller")) %>%
  mutate(Treatment = if_else(str_detect(Condition, "High"), "high", "low")) %>%
  write_csv("diff-bind-sample-sheet.csv")

samples <- dba(sampleSheet = "diff-bind-sample-sheet.csv")

myc <- dba(samples, mask = samples$masks$MYC) 
myc_counts <- dba.count(myc, summits = 250, bParallel = T)
myc_diff <- dba.contrast(myc_counts, categories = DBA_TREATMENT, minMembers = 2)
myc_res <- dba.analyze(myc_diff)
myc_res_ranges <- dba.report(myc_res)

final_anno <- annotatePeakInBatch(
  myc_res_ranges,
  AnnotationData = annoData,
  output = "overlapping",
  maxgap = 5000L
) %>%
  as.data.frame() %>%
  left_join(grch37, by = c("feature" = "ensgene"))

myc_anno_ranges <- annotatePeakInBatch(myc_res_ranges,
                                       AnnotationData = annoData,
                                       output = "overlapping",
                                       maxgap = 5000L)

genes_of_interest <- c("MYC", "EPHA2", "ETS2", "TNFAIP2", "SAMD9", "MEIS1", "IL6R", "FHL2", "C10orf10", "ISL1", "DUSP5")

sig_features <- final_anno %>%
  dplyr::filter(symbol %in% genes_of_interest)

sig_GRanges <- myc_anno_ranges[myc_anno_ranges$feature %in% sig_features$feature]
