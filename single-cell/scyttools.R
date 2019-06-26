#!/usr/bin/env Rscript

require(docopt)
require(methods)

"
Usage:
cledger.R (-h | --help | --version)
cledger.R --normalize_RSEM DIR OUT
cledger.R --normalize_featureCounts DIR OUT

Description:   This program is a command line interface to edgeR

Options:

--version                  Show the current version.
--normalize_RSEM           Perform filtering and TMM normalization on RSEM quantified RNA seq data
--normalize_featureCounts  Perform filtering and TMM normalization on featureCounts quantified RNA seq data

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
  
  if(args$`--normalize_RSEM` == T){ # still experimental
    # grab directory from command line argument
    dir_name <- args$DIR
    # recursively list all rsem output files in the directories
    file_names <- list.files(path = dir_name,
                             pattern = ".genes.results$|.genes.results.gz$",
                             full.names = T,
                             recursive = T)
    # strip out unnecessary substrings
    chr_vec <- file_names
    match_end <- matrix(nrow = length(chr_vec), ncol = length(chr_vec))
    for(str_index in 1:length(chr_vec)){
      if(str_index == length(chr_vec)) break
      string_a <- chr_vec[str_index]
      first_indx_str_b <- str_index + 1
      for(str_b_index in first_indx_str_b:length(chr_vec)){
        string_b <- chr_vec[str_b_index]
        i <- 1
        while(str_sub(string_a, end = i) == str_sub(string_b, end = i) ){
          i <- i + 1
        }
        match_end[str_index, str_b_index] <- i
      }
    }
    # create vector or new column names
    new_col_names <- file_names %>%
      map(str_sub,
          start = min(match_end, na.rm = T)) %>%
      unlist() %>%
      str_remove_all("\\.genes\\.results|_star_rsem")
    # perform tximport
    txi <- tximport(file_names, type = "rsem", txIn = FALSE, txOut = FALSE)
    txi$length[txi$length == 0] <- 1
    # perform offset calulation
    cts <- txi$counts
    normMat <- txi$length
    normMat <- normMat/exp(rowMeans(log(normMat)))
    
    o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))
    # construct edger object
    y <- DGEList(cts)
    y <- scaleOffset(y, t(t(log(normMat)) + o))
    # filtering
    keep <- filterByExpr(y)
    y <- y[keep, ]
    # normalize
    y <- calcNormFactors(y)
    file_out <- args$OUT
    # write out transformed CPM values
    cpm(y, log = T) %>%
      as.data.frame() %>%
      rownames_to_column("gene_id") %>%
      setNames(c("gene_id", new_col_names)) %>%
      write_csv(file_out)
  }else if (args$`--normalize_featureCounts` == T){
    dir_name <- args$DIR
    # recursively list all rsem output files in the directories
    file_names <- list.files(path = dir_name,
                             pattern = ".txt$",
                             full.names = T,
                             recursive = F)
    
    cts <- map(as.character(file_names), read_tsv, skip = 1) %>%
      purrr::reduce(left_join)
    
    # construct edger object
    y <- DGEList(counts = cts %>%
                   dplyr::select(ends_with("sortedByCoord.out.bam")),
                 genes = cts$Geneid,
                 remove.zeros = T)
    
    # filtering
    keep <- filterByExpr(y)
    y <- y[keep, ]
    # normalize
    y <- calcNormFactors(y)
    file_out <- args$OUT
    # write out transformed CPM values
    cpm(y, log = T) %>%
      as.data.frame() %>%
      bind_cols(y$genes, .) %>%
      write_csv(file_out)
  }else{
    cat(paste(c("\nERROR:","\nCommand not found"), collapse = "\n"), "\n")  
  }
}
