#' ---
#' title: "Playing with spin"
#' author: "Brian Capaldo"
#' date: "3/05/2019"
#' output: github_document
#' ---

#' # load libraries

#+ r setup, include=FALSE

library(edgeR)
library(tidyverse)
library(clusterProfiler)
library(annotables)
library(broom)
library(scran)
library(scater)
library(SC3)
library(BiocNeighbors)
library(BiocParallel)
library(DropletUtils)

#' # analysis
#' 
#' Running vignette developed by Stephanie Hicks currently available [here](https://bioconductor.github.io/OrchestratingSingleCellAnalysis/workflow-integrating-datasets.html)

#' ## building out filitering criteria
#' 
#' Mitochondria genes, cell cycle genes, and ribosomal genes are the most commonly discussed sets of genes that can impact single cell analysis. To address this, we create three annotated tables for each of these sets. Ribosomal genes are defined as those with symbols starting with RPL or RPS and are protein coding. Mitochondria genes are defined as those with symbols starting with MT- and are protein coding. Cell Cycle genes have been previously defined by INSERT REFERENCE. We also take the opportunity to load in the MSigDB, as it has many genesets that will serve as sanity checks throughout the analysis. 

# ribosomal genes
rb_genes <- grch38 %>%
  dplyr::filter(str_detect(symbol, "^RPL|^RPS"),
                biotype == "protein_coding") 

# mitchondrial genes
mt_genes <- grch38 %>%
  dplyr::filter(str_detect(symbol, "^MT-"),
                biotype == "protein_coding") 

# load msigdb
msig_db <- read.gmt("/Volumes/Group05/CCBB/cBioPortal-Projects/msigdb.v6.2.entrez.gmt")


#' # references