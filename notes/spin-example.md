Playing with spin
================
Brian Capaldo
3/05/2019

load libraries
==============

analysis
========

Running vignette developed by Stephanie Hicks currently available [here](https://bioconductor.github.io/OrchestratingSingleCellAnalysis/workflow-integrating-datasets.html) \#\# building out filitering criteria

Mitochondria genes, cell cycle genes, and ribosomal genes are the most commonly discussed sets of genes that can impact single cell analysis. To address this, we create three annotated tables for each of these sets. Ribosomal genes are defined as those with symbols starting with RPL or RPS and are protein coding. Mitochondria genes are defined as those with symbols starting with MT- and are protein coding. Cell Cycle genes have been previously defined by INSERT REFERENCE. We also take the opportunity to load in the MSigDB, as it has many genesets that will serve as sanity checks throughout the analysis.

``` r
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
```

    ## Loading required package: GSEABase

    ## Loading required package: annotate

    ## Loading required package: AnnotationDbi

    ## 
    ## Attaching package: 'AnnotationDbi'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select

    ## Loading required package: XML

    ## Warning: package 'XML' was built under R version 3.5.2

    ## 
    ## Attaching package: 'annotate'

    ## The following object is masked from 'package:BiocNeighbors':
    ## 
    ##     findNeighbors

    ## Loading required package: graph

    ## 
    ## Attaching package: 'graph'

    ## The following object is masked from 'package:XML':
    ## 
    ##     addNode

    ## The following object is masked from 'package:stringr':
    ## 
    ##     boundary

    ## Warning in getGmt(con = gmtfile): 1 record(s) contain duplicate ids:
    ## QUINTENS_EMBRYONIC_BRAIN_RESPONSE_TO_IR

references
==========
