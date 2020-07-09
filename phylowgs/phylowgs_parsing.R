#!/usr/bin/env Rscript

require(docopt)
require(methods)

"
Usage:
phylowgs_parsing.R (-h | --help | --version)
phylowgs_parsing.R DIR OUT

Description:   This script parses phyloWGS output

Options:
--version       Show the current version.

Arguments:
DIR    Provide directory where phyloWGS output is located

" -> doc


args <- docopt(doc)

dir <- args$DIR

library(tidyverse)
library(jsonlite)
library(igraph)
library(tidygraph)
library(ggraph)
library(patchwork)

setwd(dir)

cnv_tar <- untar("cnv_data.tar",
                 list = T)

ssm_files <- untar("ssm_data.tar",
                   list = T)

summary_files <- list.files(pattern = ".summ.json",
                            recursive = T,
                            full.names = T)

best_trees_list <- lapply(summary_files, function(summary_file){
  cat(summary_file, "\n")
  # load in summary data from phylowgs as list
  summary_data <- jsonlite::fromJSON(summary_file)
  
  # identify tree with the maximum log likelihood (llh)
  tree_id <- lapply(summary_data$trees, function(tree){
    return(data.frame(llh = tree$llh))
  }) %>%
    bind_rows(.id = "tree_id") %>%  
    filter(llh == max(llh))
  
  # load tables with ssm/cna information
  tree_alt_files_df <- tree_id %>% 
    mutate(summary_file = c(summary_file),
           tree_file = paste(tree_id, "json", sep = "."),
           zip_archive = str_replace(summary_file, "summ.json", "mutass.zip"))
  
  unzip(tree_alt_files_df$zip_archive[1], tree_alt_files_df$tree_file[1], exdir = "unzipped_mutass_trees")
  tree <- fromJSON(paste0("unzipped_mutass_trees/", tree_alt_files_df$tree_file[1]))
  # make data.frame of alts
  alts_df <- lapply(tree$mut_assignments, function(mut_assignment){
    return(data.frame(ssm_cnv = unlist(c(mut_assignment$ssms, mut_assignment$cnvs))))
  }) %>% bind_rows(.id = "node_id") %>% 
    mutate(tree_id = rep(tree_id$tree_id, nrow(.)))
  
  # subset out the tree with the maximum (llh)
  tree_info <- summary_data$trees[[tree_id$tree_id]]
  
  # extract structure into data.frame
  structure_df <- tree_info$structure %>%
    lapply(function(parent){
      return(data.frame(child = parent))
    }) %>% 
    bind_rows(.id = "parent") %>% 
    mutate(tree_id = rep(tree_id$tree_id[1], nrow(.)))
  
  # extract population information into data.frame
  populations_df <- tree_info$populations %>%
    lapply(function(node){
      node_df <- tibble(num_ssms = node$num_ssms,
                        num_cnvs = node$num_cnvs,
                        cellular_prevalence = node$cellular_prevalence) %>%
        mutate(sample_index = paste0("sample_", c(1:nrow(.)))) %>% 
        nest(cellular_prevalence = c(cellular_prevalence, sample_index))
    }) %>% 
    bind_rows(.id = "node_id") %>% 
    mutate(tree_id = rep(tree_id$tree_id[1], nrow(.)))
  
  # prune trees
  # each node (or population) must have at least 5 ssms or cnas
  #  nodes without siblings are merged into parent nodes
  #  nodes with siblings are merged into their siblings
  #  nodes are eliminated if they are the only child of the normal node (node_id of parent == 0)
  #   AND the children of the node become the direct descendents of the normal
  
  # create "pruned" version of populations_df, iterate, and update alts_df and structure_df as needed
  pruned_populations_df <- populations_df 
  
  pruned_alts_df <- alts_df
  
  pruned_structure_df <- structure_df
  
  # ensure that the normal node has multiple children
  
  if(table(pruned_structure_df$parent)["0"] == 1){
    
    big_brother <- "0"
    node_to_collapse <- "1"
    
    # update propulations_df
    node_updated_pops_df <- pruned_populations_df %>% 
      mutate(node_id = as.character(node_id),
             node_id = if_else(node_id == node_to_collapse, as.character(big_brother), node_id)) %>%
      filter(node_id == big_brother) %>% 
      unnest(cellular_prevalence) %>% 
      group_by(node_id, sample_index) %>% 
      summarize(num_ssms = sum(num_ssms),
                num_cnvs = sum(num_cnvs),
                cellular_prevalence = cellular_prevalence[2],
                tree_id = max(tree_id)) %>% 
      ungroup() 
    
    pruned_populations_df <- pruned_populations_df %>% 
      mutate(node_id = as.character(node_id),
             node_id = if_else(node_id == node_to_collapse, as.character(big_brother), node_id)) %>%
      filter(node_id != big_brother) %>%
      bind_rows(node_updated_pops_df %>% 
                  nest(cellular_prevalence = c(cellular_prevalence, sample_index))) %>% 
      mutate(node_id = factor(node_id, levels = c(0:100))) %>% 
      arrange(node_id) %>% 
      mutate(node_id = as.character(node_id))
    
    # update alts_df
    pruned_alts_df <- pruned_alts_df %>% 
      mutate(node_id = if_else(node_id == node_to_collapse, as.character(big_brother), node_id))
    # update structure_df
    # update structure_df
    pruned_structure_df <- pruned_structure_df %>% 
      mutate(parent = if_else(parent == node_to_collapse, as.character(big_brother), parent),
             child = if_else(child == node_to_collapse, as.character(big_brother), as.character(child))) 
    
    if(nrow(pruned_structure_df) > 1){
      pruned_structure_df <- pruned_structure_df %>% 
        filter(parent != child) %>% 
        distinct()
    }
  }
  
  # identify most terminal node to collapse
  node_to_collapse <- (pruned_populations_df %>% 
                         filter(num_cnvs < 5,
                                num_ssms < 5,
                                node_id != 0) %>% 
                         slice(nrow(.)))$node_id
  # if a node needs to be collapse, proceed to collapse, otherwise do nothing
  while(length(node_to_collapse) > 0){
    # identify parent
    parent_node <- pruned_structure_df[pruned_structure_df$child == node_to_collapse,]$parent 
    
    # identify siblings
    sibling_nodes <- pruned_structure_df[pruned_structure_df$parent == parent_node,]$child 
    sibling_nodes <- sibling_nodes[sibling_nodes != node_to_collapse]
    
    # check for the existance of siblings, if true, find big brother, if not, assign parent node as big brother
    if(length(sibling_nodes) > 0){
      # assign big brother, or the node to which sibling nodes will be collapsed to (this selects the node with the most alterations)
      big_brother <- (pruned_populations_df %>% 
                        filter(node_id %in% sibling_nodes) %>% 
                        arrange(num_ssms) %>% 
                        slice(nrow(.)))$node_id
      
      # update propulations_df
      node_updated_pops_df <- pruned_populations_df %>% 
        mutate(node_id = as.character(node_id),
               node_id = if_else(node_id == node_to_collapse, as.character(big_brother), node_id)) %>%
        filter(node_id == big_brother) %>% 
        unnest(cellular_prevalence) %>% 
        group_by(node_id, sample_index) %>% 
        summarize(num_ssms = sum(num_ssms),
                  num_cnvs = sum(num_cnvs),
                  cellular_prevalence = sum(cellular_prevalence),
                  tree_id = max(tree_id)) %>% 
        ungroup() 
      
      pruned_populations_df <- pruned_populations_df %>% 
        mutate(node_id = as.character(node_id),
               node_id = if_else(node_id == node_to_collapse, as.character(big_brother), node_id)) %>%
        filter(node_id != big_brother) %>%
        bind_rows(node_updated_pops_df %>% 
                    nest(cellular_prevalence = c(cellular_prevalence, sample_index))) %>% 
        mutate(node_id = factor(node_id, levels = c(0:100))) %>% 
        arrange(node_id) %>% 
        mutate(node_id = as.character(node_id))
    }else{
      big_brother <- parent_node
      
      # update propulations_df
      node_updated_pops_df <- pruned_populations_df %>% 
        mutate(node_id = as.character(node_id),
               node_id = if_else(node_id == node_to_collapse, as.character(big_brother), node_id)) %>%
        filter(node_id == big_brother) %>% 
        unnest(cellular_prevalence) %>% 
        group_by(node_id, sample_index) %>% 
        summarize(num_ssms = sum(num_ssms),
                  num_cnvs = sum(num_cnvs),
                  cellular_prevalence = cellular_prevalence[1],
                  tree_id = max(tree_id)) %>% 
        ungroup() 
      
      pruned_populations_df <- pruned_populations_df %>% 
        mutate(node_id = as.character(node_id),
               node_id = if_else(node_id == node_to_collapse, as.character(big_brother), node_id)) %>%
        filter(node_id != big_brother) %>%
        bind_rows(node_updated_pops_df %>% 
                    nest(cellular_prevalence = c(cellular_prevalence, sample_index))) %>% 
        mutate(node_id = factor(node_id, levels = c(0:100))) %>% 
        arrange(node_id) %>% 
        mutate(node_id = as.character(node_id))
    }
    # update alts_df
    pruned_alts_df <- pruned_alts_df %>% 
      mutate(node_id = if_else(node_id == node_to_collapse, as.character(big_brother), node_id))
    # update structure_df
    pruned_structure_df <- pruned_structure_df %>% 
      mutate(parent = if_else(parent == node_to_collapse, as.character(big_brother), parent),
             child = if_else(child == node_to_collapse, as.character(big_brother), as.character(child))) 
    
    if(nrow(pruned_structure_df) > 1){
      pruned_structure_df <- pruned_structure_df %>% 
        filter(parent != child) %>% 
        distinct()
    }
    # update node_to_collapse
    node_to_collapse <- (pruned_populations_df %>% 
                           filter(num_cnvs < 5,
                                  num_ssms < 5,
                                  node_id != 0) %>% 
                           slice(nrow(.)))$node_id
  }
  
  tree_layout <- graph_from_data_frame(structure_df) %>% 
    create_layout("tree") %>% 
    mutate(tree_level = factor(y, levels = c(0:max(y))))
  cancer_cell_fraction <- populations_df %>% 
    unnest(cols = c(cellular_prevalence)) %>% 
    mutate(ccf = rep(0, nrow(.)))
  for(tree_level in levels(tree_layout$tree_level)){
    for(node in tree_layout$name[tree_layout$tree_level == tree_level]){
      children <- structure_df$child[structure_df$parent == node]
      if(length(children) == 0){
        cancer_cell_fraction <- cancer_cell_fraction %>% 
          mutate(ccf = if_else(node_id == node, cellular_prevalence, ccf))
      }else{
        children_ccf <- cancer_cell_fraction %>% 
          filter(node_id %in% children) %>% 
          group_by(sample_index) %>% 
          summarise(sub_ccf = sum(cellular_prevalence)) %>% 
          ungroup() %>% 
          mutate(node_id = rep(as.character(node), nrow(.)))
        
        cancer_cell_fraction <- cancer_cell_fraction %>% 
          left_join(children_ccf) %>% 
          mutate(ccf = if_else(node_id == node, cellular_prevalence - sub_ccf, ccf)) %>% 
          dplyr::select(-sub_ccf)
      }
    }
  }
  cancer_cell_fraction$ccf[cancer_cell_fraction$node_id == 0] <- NA
  
  tree_layout <- graph_from_data_frame(pruned_structure_df) %>% 
    create_layout("tree") %>% 
    mutate(tree_level = factor(y, levels = c(0:max(y))))
  pruned_cancer_cell_fraction <- pruned_populations_df %>% 
    unnest(cols = c(cellular_prevalence)) %>% 
    mutate(ccf = rep(0, nrow(.)))
  for(tree_level in levels(tree_layout$tree_level)){
    for(node in tree_layout$name[tree_layout$tree_level == tree_level]){
      children <- pruned_structure_df$child[pruned_structure_df$parent == node]
      if(length(children) == 0){
        pruned_cancer_cell_fraction <- pruned_cancer_cell_fraction %>% 
          mutate(ccf = if_else(node_id == node, cellular_prevalence, ccf))
      }else{
        children_ccf <- pruned_cancer_cell_fraction %>% 
          filter(node_id %in% children) %>% 
          group_by(sample_index) %>% 
          summarise(sub_ccf = sum(cellular_prevalence)) %>% 
          ungroup() %>% 
          mutate(node_id = rep(as.character(node), nrow(.)))
        
        pruned_cancer_cell_fraction <- pruned_cancer_cell_fraction %>% 
          left_join(children_ccf) %>% 
          mutate(ccf = if_else(node_id == node, cellular_prevalence - sub_ccf, ccf)) %>% 
          dplyr::select(-sub_ccf)
      }
    }
  }
  if(any(unique(pruned_cancer_cell_fraction$node_id) == "1")){
    pruned_cancer_cell_fraction$ccf[pruned_cancer_cell_fraction$node_id == 0] <- NA
  }
  
  shannon_diversity_index <- cancer_cell_fraction %>% 
    dplyr::filter(!is.na(ccf)) %>% 
    mutate(to_be_summed = ccf*log(ccf)) %>% 
    group_by(tree_id, sample_index) %>% 
    summarise(div_index = -sum(to_be_summed))
  
  pruned_shannon_diversity_index <- pruned_cancer_cell_fraction %>% 
    dplyr::filter(!is.na(ccf)) %>% 
    mutate(to_be_summed = ccf*log(ccf)) %>% 
    group_by(tree_id, sample_index) %>% 
    summarise(div_index = -sum(to_be_summed))
  
  tree_id <- tree_id %>%
    left_join(structure_df %>%
                nest(structure = c(parent, child)),
              by = "tree_id")  %>%
    left_join(pruned_structure_df %>%
                nest(pruned_structure = c(parent, child)),
              by = "tree_id")  %>%
    left_join(populations_df %>%
                nest(populations = c(node_id:cellular_prevalence)),
              by = "tree_id")  %>% 
    left_join(pruned_populations_df %>%
                nest(pruned_populations = c(node_id:cellular_prevalence)),
              by = "tree_id")  %>% 
    left_join(alts_df %>% 
                nest(alts = c(node_id:ssm_cnv)),
              by = "tree_id")  %>% 
    left_join(pruned_alts_df %>% 
                nest(pruned_alts = c(node_id:ssm_cnv)),
              by = "tree_id") %>% 
    left_join(cancer_cell_fraction %>% 
                dplyr::select(tree_id, everything()) %>% 
                nest(cancer_cell_fraction = c(node_id:ccf)),
              by = "tree_id") %>% 
    left_join(pruned_cancer_cell_fraction %>% 
                dplyr::select(tree_id, everything()) %>% 
                nest(pruned_cancer_cell_fraction = c(node_id:ccf)),
              by = "tree_id") %>% 
    left_join(shannon_diversity_index %>% 
                nest(shannon_diversity_index = c(sample_index:div_index)),
              by = "tree_id") %>% 
    left_join(pruned_shannon_diversity_index %>% 
                nest(pruned_shannon_diversity_index = c(sample_index:div_index)),
              by = "tree_id")
  
  return(tree_id)  
})

names(best_trees_list) <- summary_files

best_trees <- best_trees_list %>% 
  bind_rows(.id = "summary_file") %>% 
  mutate(sample_id = str_remove(summary_file, ".summ.json$"))

dir.create("alt_tables", showWarnings = F)
dir.create("phylogency-plots", showWarnings = F)

for(sample in unique(best_trees$sample_id)){
  cat(sample, "\n")
  plot_trees <- best_trees %>% 
    filter(sample_id == sample)
    plot_tree <- plot_trees 
    p1 <- plot_tree$structure[[1]] %>% 
      graph_from_data_frame() %>% 
      as_tbl_graph() %>% 
      activate(nodes) %>% 
      left_join(plot_tree$cancer_cell_fraction[[1]] %>% 
                  group_by(node_id) %>% 
                  summarise(ccf = mean(ccf, na.rm = T)),
                by = c("name" = "node_id")) %>% 
      create_layout("tree") %>% 
      ggraph() +
      geom_node_text(aes(label = name), repel = T) + 
      geom_node_point(aes(size = ccf)) +
      geom_edge_link() +
      theme_void() +
      theme(legend.position = "none")
    
    p2 <- plot_tree$pruned_structure[[1]] %>% 
      graph_from_data_frame() %>% 
      as_tbl_graph() %>% 
      activate(nodes) %>% 
      left_join(plot_tree$pruned_cancer_cell_fraction[[1]] %>% 
                  group_by(node_id) %>% 
                  summarise(ccf = mean(ccf, na.rm = T)),
                by = c("name" = "node_id")) %>% 
      create_layout("tree") %>% 
      ggraph() +
      geom_node_text(aes(label = name), repel = T) + 
      geom_node_point(aes(size = ccf)) +
      theme_void() +
      theme(legend.position = "none")
    
    if(nrow(plot_tree$pruned_structure[[1]]) > 1){
      p2 <- p2 +
        geom_edge_link()
      }
    
    p1+p2
    
    plot_file_out <- paste0("phylogency-plots/", sample, "_", ".svg")
    
    ggsave(plot_file_out, height = 9, width = 16)
    
    untar("ssm_data.tar", files = paste(sample, "ssm_data.txt", sep = "/"))
    
    plot_tree$alts[[1]] %>%
      filter(str_detect(ssm_cnv, "^s")) %>%
      left_join(read_tsv(paste(sample, "ssm_data.txt", sep = "/")),
                by = c("ssm_cnv" = "id")) %>% 
      write_csv(paste0("alt_tables/", sample, "_ssm_table.csv"))
    
    plot_tree$pruned_alts[[1]] %>%
      filter(str_detect(ssm_cnv, "^s")) %>%
      left_join(read_tsv(paste(sample, "ssm_data.txt", sep = "/")),
                by = c("ssm_cnv" = "id")) %>% 
      write_csv(paste0("alt_tables/", sample, "_pruned_ssm_table.csv"))
    
    untar("cnv_data.tar", files = paste(sample, "cnv_data.txt", sep = "/"))
    
    plot_tree$alts[[1]] %>%
      filter(str_detect(ssm_cnv, "^c")) %>%
      left_join(read_tsv(paste(sample, "cnv_data.txt", sep = "/")),
                by = c("ssm_cnv" = "cnv")) %>% 
      write_csv(paste0("alt_tables/", sample, "_cnv_table.csv"))
    
    plot_tree$pruned_alts[[1]] %>%
      filter(str_detect(ssm_cnv, "^c")) %>%
      left_join(read_tsv(paste(sample, "cnv_data.txt", sep = "/")),
                by = c("ssm_cnv" = "cnv")) %>% 
      write_csv(paste0("alt_tables/", sample, "_pruned_cnv_table.csv"))
}

best_trees %>% 
  unnest(pruned_shannon_diversity_index) %>% 
  group_by(sample_id) %>% 
  summarize(shannon_index = mean(div_index)) %>% 
  write_csv("shannon-index-table.csv")
