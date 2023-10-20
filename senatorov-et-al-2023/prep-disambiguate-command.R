library(tidyverse)

human_alignments <- list.files("/data/capaldobj/incoming-nih-dme/CS035088-chip-seq-ilya-hnf1a/chipseq-results/STAR/",
full.names = T, 
pattern = ".bam$")

mouse_alignments <- list.files("/data/capaldobj/CS026803-rna-seq-yin-ivy-bone-marrow-metastatis/grcm38-rna-seq-results/STAR/",
full.names = T, 
pattern = ".bam$")

human_df <- data.frame(bam_file_human = as.character(human_alignments),
                       sample_id = basename(human_alignments))
mouse_df <- data.frame(bam_file_mouse = as.character(mouse_alignments),
                       sample_id = basename(mouse_alignments))

human_df %>% 
  left_join(mouse_df) %>% 
  mutate(sample_id = str_remove(sample_id, "_S[0-9]*_R1_001Aligned.sortedByCoord.out.bam"),
         command = paste0("ngs_disambiguate -s ", sample_id, " -o ngs_disambiguated -a star ",
                          bam_file_human,
                          " ",
                          bam_file_mouse)) %>% 
  select(command) %>% 
  write_tsv("/data/capaldobj/CS026803-rna-seq-yin-ivy-bone-marrow-metastatis/ngs_disambiguate.swarm", col_names = F)

data.frame(bam_file = list.files("/data/capaldobj/CS026803-rna-seq-yin-ivy-bone-marrow-metastatis/ngs_disambiguated/",
                                  pattern = ".disambiguatedSpeciesA.bam",
                                  full.names = T),
           sample_id = basename(list.files("/data/capaldobj/CS026803-rna-seq-yin-ivy-bone-marrow-metastatis/ngs_disambiguated/",
                                  pattern = ".disambiguatedSpeciesA.bam",
                                  full.names = T))) %>% 
   mutate(sample_id = str_remove(sample_id, ".disambiguatedSpeciesA.bam"),
         out_file = paste0("/data/capaldobj/CS026803-rna-seq-yin-ivy-bone-marrow-metastatis/featureCounts/", sample_id, "_grch37_counts.txt"),
         command = paste0("featureCounts -T 4 -t exon -g gene_id -a /fdb/igenomes/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf -o ",
                          out_file,
                          " ",
                          bam_file))%>% 
  select(command) %>% 
  write_tsv("/data/capaldobj/CS026803-rna-seq-yin-ivy-bone-marrow-metastatis/featureCounts-human.swarm", col_names = F)

data.frame(bam_file = list.files("/data/capaldobj/CS026803-rna-seq-yin-ivy-bone-marrow-metastatis/ngs_disambiguated/",
                                  pattern = ".disambiguatedSpeciesB.bam",
                                  full.names = T),
           sample_id = basename(list.files("/data/capaldobj/CS026803-rna-seq-yin-ivy-bone-marrow-metastatis/ngs_disambiguated/",
                                  pattern = ".disambiguatedSpeciesB.bam",
                                  full.names = T))) %>% 
   mutate(sample_id = str_remove(sample_id, ".disambiguatedSpeciesB.bam"),
         out_file = paste0("/data/capaldobj/CS026803-rna-seq-yin-ivy-bone-marrow-metastatis/featureCounts/", sample_id, "_grcm38_counts.txt"),
         command = paste0("featureCounts -T 4 -t exon -g gene_id -a /fdb/igenomes/Mus_musculus/Ensembl/GRCm38/Annotation/Genes/genes.gtf -o ",
                          out_file,
                          " ",
                          bam_file))%>% 
  select(command) %>% 
  write_tsv("/data/capaldobj/CS026803-rna-seq-yin-ivy-bone-marrow-metastatis/featureCounts-mouse.swarm", col_names = F)