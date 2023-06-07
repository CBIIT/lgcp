library(tidyverse)

dir.create("awesome-prostate-signatures")

read_csv("~/lgcp/rnaseq/analysis-scripts/ar_reference_beltran_cpm_grch37.csv") %>% 
  write_csv("awesome-prostate-signatures/ar_reference_beltran_cpm_grch37.csv")

read_csv("~/lgcp/rnaseq/analysis-scripts/neuro_reference_beltran_cpm_grch37.csv") %>% 
  write_csv("awesome-prostate-signatures/neuro_reference_beltran_cpm_grch37.csv")

read_tsv("/Volumes/Group05/CCBB/CS025602_PDO_PDX/metadata-from-lab/PRAD.norm_Beltran_LUAD.norm_SCLC_LUAD.subset_rsem_genes_upper_norm_counts_coding_log2_rmprolifROC_prcomp_loadings_VARIMAX.txt") %>% 
  write_tsv("awesome-prostate-signatures/PRAD.norm_Beltran_LUAD.norm_SCLC_LUAD.subset_rsem_genes_upper_norm_counts_coding_log2_rmprolifROC_prcomp_loadings_VARIMAX.txt")

read_tsv("/Volumes/Group05/CCBB/CS025602_PDO_PDX/metadata-from-lab/PRAD.norm_Beltran_LUAD.norm_SCLC_LUAD.subset_rsem_genes_upper_norm_counts_coding_log2_prolifgenesROC_prcomp_loadings.txt") %>% 
  write_tsv("awesome-prostate-signatures/PRAD.norm_Beltran_LUAD.norm_SCLC_LUAD.subset_rsem_genes_upper_norm_counts_coding_log2_prolifgenesROC_prcomp_loadings.txt")

read_csv("/Volumes/Group05/CCBB/Single-Cell-Bioinformatics-2019-October-03/signature-table-2023-May-01.csv") %>% 
  write_csv("awesome-prostate-signatures/signature-table-2023-May-01.csv")
