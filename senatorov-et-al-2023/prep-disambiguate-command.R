library(tidyverse)

human_alignments <- list.files("/data/capaldobj/incoming-nih-dme/CS035088-chip-seq-ilya-hnf1a/chipseq-grch37-results/bwa/mergedLibrary/",
full.names = T,
pattern = ".bam$")

mouse_alignments <- list.files("/data/capaldobj/incoming-nih-dme/CS035088-chip-seq-ilya-hnf1a/chipseq-results/bwa/mergedLibrary/",
full.names = T,
pattern = ".bam$")

human_df <- data.frame(bam_file_human = as.character(human_alignments),
                       sample_id = basename(human_alignments))
mouse_df <- data.frame(bam_file_mouse = as.character(mouse_alignments),
                       sample_id = basename(mouse_alignments))

human_df %>%
  left_join(mouse_df) %>%
  mutate(sample_id = str_remove(sample_id,
    ".mLb.clN.sorted.bam"),
         command = paste0("ngs_disambiguate -s ",
         sample_id, " -o ngs_disambiguated -a bwa ",
                          bam_file_human,
                          " ",
                          bam_file_mouse)) %>%
  select(command) %>%
  write_tsv("/data/capaldobj/lgcp/senatorov-et-al-2023/ngs_disambiguate.swarm", col_names = F)

# swarm -f /data/capaldobj/lgcp/senatorov-et-al-2023/ngs_disambiguate.swarm -t 8 -g 64 --partition ccr