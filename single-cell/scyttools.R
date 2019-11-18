#!/usr/bin/env Rscript

require(docopt)
require(methods)

"
Usage:
scyttools.R (-h | --help | --version)
scyttools.R --quality_control DIR OUT
scyttools.R --dimensionality_reduction RDS OUT
scyttools.R --cluster_cells RDS OUT
scyttools.R --trajectory_inference RDS OUT
scyttools.R --geneset_scoring RDS OUT
scyttools.R --differential_expression RDS OUT

Description:   This program is a command line interface to running automated single cell RNA sequencing analysis

Options:

--version                     Show the current version.
--quality_control             Perform quality control and remove cells leveraging the Scater and Scran packages
--dimensionality_reduction    Perform dimensionality reduction using GLM PCA
--cluster_cells               Cluster cells using GLM PCA and self-organizing maps
--trajectory_inference        Perform trajectory inference using monocle version 2
--geneset_scoring             Score cells and clusters for geneset activity
--differential_expression     Perform differential expression analysis of clusters

Arguments:

DIR   Directory that contains cellranger gene barcode matrix output
OUT   Provide output file
RDS   RDS file that contains saved SingleCellExperiment object


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
