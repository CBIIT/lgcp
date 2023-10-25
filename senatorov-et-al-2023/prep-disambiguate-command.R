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
  write_tsv("/data/capaldobj/lgcp/senatorov-et-al-2023/ngs_disambiguate.swarm",
  col_names = F)

# swarm -f /data/capaldobj/lgcp/senatorov-et-al-2023/ngs_disambiguate.swarm -t 8 -g 64 --partition ccr

# convert 170.2 disambiguated files and normal 23.1 alignments back to fastq

# grab disambiguated 170 bams
disambiguated_bams <- list.files(
  "/data/capaldobj/incoming-nih-dme/CS035088-chip-seq-ilya-hnf1a/ngs_disambiguated/",
  full.names = T)

disambiguated_bams <- disambiguated_bams[str_detect(disambiguated_bams, "170")]
disambiguated_bams <- disambiguated_bams[
  str_detect(disambiguated_bams,
  ".disambiguatedSpeciesA.bam")]

# combine with 23.1 bams
data.frame(disambiguated_bams = c(human_alignments[
  !str_detect(human_alignments, "170")],
  disambiguated_bams)) %>%
  mutate(command = paste0("samtools sort -n ",
    disambiguated_bams,
    " -o ",
    str_replace(basename(disambiguated_bams),
      "mLb.clN.sorted|disambiguatedSpeciesA",
      "sorted"),
    "; ",
    "samtools fastq -@ 8 ",
    str_replace(basename(disambiguated_bams),
      "mLb.clN.sorted|disambiguatedSpeciesA",
      "sorted"),
    " -1 ",
    str_replace(basename(disambiguated_bams),
      "mLb.clN.sorted.bam|disambiguatedSpeciesA.bam",
      "_R1.fastq.gz"),
    " -2 ",
    str_replace(basename(disambiguated_bams),
      "mLb.clN.sorted.bam|disambiguatedSpeciesA.bam",
      "_R2.fastq.gz"),
    " -0 /dev/null -s /dev/null -n")
  ) %>%
  select(command) %>%
  write_tsv("/data/capaldobj/lgcp/senatorov-et-al-2023/bam-to-fastq.swarm",
  col_names = F)

# generate new sample sheet
data.frame(full_file_name = list.files("/data/capaldobj/incoming-nih-dme/CS035088-chip-seq-ilya-hnf1a/bam2fastq-results",
  full.names = T,
  pattern = ".fastq.gz")) %>%
  mutate(sample_id = basename(full_file_name) %>%
    str_remove("\\._R[1-2].fastq.gz"),
            fastq_id = if_else(str_detect(full_file_name, "R2"),
                "fastq_2",
                "fastq_1")) %>%
    pivot_wider(names_from = fastq_id,
        values_from = full_file_name) %>%
        unnest(cols = c(fastq_1, fastq_2)) %>%
    mutate(antibody = case_when(str_detect(sample_id, "HNF1A") ~ "HNF1A",
      str_detect(sample_id, "H3K27Ac") ~ "H3K27Ac",
      T ~ ""),
      control = if_else(str_detect(sample_id, "input"),
        "",
        str_replace(sample_id,
        "_HNF1A|_H3K27Ac",
        "_TE_input"))) %>%
      arrange(desc(antibody)) %>%
      setNames(c("sample","fastq_1","fastq_2","antibody","control")) %>%
    write_csv("senatorov-et-al-2023/design-dis-hnf1a-chip.csv")
