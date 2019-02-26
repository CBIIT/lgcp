#!/usr/bin/env Rscript

require(docopt)
require(methods)

"
Usage:
cledger.R (-h | --help | --version)
cledger.R --normalize RSEM OUT

Description:   This program is a command line interface to edgeR

Options:

--version   Show the current version.
--normalize Perform filtering and TMM normalization on quantified RNA seq data

Arguments:

RSEM  Provide directory with RSEM .genes.results file(s)
OUT   Provide output file

" -> doc

args <- docopt(doc)

if(args$`--version` == T){ # returns version if version is requested
  cat("\nVersion is 0.1\n")
}else{ # Analysis begins
  
  library(edgeR)
  library(tidyverse)
  library(tximport)
  
  if(args$`--normalize`){
    dir_name <- args$RSEM
    file_names <- list.files(path = dir_name,
                             pattern = ".genes.results$|.genes.results.gz$",
                             full.names = T,
                             recursive = T)
    txi <- tximport(file_names, type = "rsem", txIn = FALSE, txOut = FALSE)
    txi$length[txi$length == 0] <- 1
    
    cts <- txi$counts
    normMat <- txi$length
    normMat <- normMat/exp(rowMeans(log(normMat)))
    
    o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))
    y <- DGEList(cts)
    y <- scaleOffset(y, t(t(log(normMat)) + o))
    # filtering
    keep <- filterByExpr(y)
    y <- y[keep, ]
    file_out <- args$OUT
    cpm(y, log = T) %>%
      as.data.frame() %>%
      rownames_to_column("gene_id") %>%
      write_csv(file_out)
  }else{
    cat(paste(c("\nERROR:","\nCommand not found"), collapse = "\n"), "\n")  
  }
}
