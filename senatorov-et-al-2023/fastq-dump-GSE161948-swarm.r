library(tidyverse)

acc_to_sample_df <- read_csv(
    "senatorov-et-al-2023/GSE161948-accession-table.txt")
# all reads are paired end, 150 bp

read_csv(
    "senatorov-et-al-2023/GSE161948-PRJNA679976-SraRunTable.csv") %>%
    mutate(command = paste(
        "fastq-dump --gzip --split-files -O /lscratch/$SLURM_JOBID",
            paste0(Run, ";"),
            "cp -R /lscratch/$SLURM_JOBID/*.fastq.gz /data/LGCP/freedman-chip/sra-files/",
            sep = " ")) %>%
    select(command) %>%
    write.table("/data/LGCP/freedman-chip/sra-files/fastq-dump-GSE161948.swarm",
        quote = F,
        row.names = F,
        col.names = F)
# swarm -f sra-files/fastq-dump.swarm --module sratoolkit --gres=lscratch:100 -g 4  -t 6 --partition=ccr

# sample,fastq_1,fastq_2,antibody,control
# WT_BCATENIN_IP_REP1,BLA203A1_S27_L006_R1_001.fastq.gz,,BCATENIN,WT_INPUT
# WT_BCATENIN_IP_REP2,BLA203A25_S16_L001_R1_001.fastq.gz,,BCATENIN,WT_INPUT
# WT_BCATENIN_IP_REP2,BLA203A25_S16_L002_R1_001.fastq.gz,,BCATENIN,WT_INPUT
# WT_BCATENIN_IP_REP2,BLA203A25_S16_L003_R1_001.fastq.gz,,BCATENIN,WT_INPUT
# WT_BCATENIN_IP_REP3,BLA203A49_S40_L001_R1_001.fastq.gz,,BCATENIN,WT_INPUT
# WT_INPUT_REP1,BLA203A6_S32_L006_R1_001.fastq.gz,,,
# WT_INPUT_REP2,BLA203A30_S21_L001_R1_001.fastq.gz,,,
# WT_INPUT_REP2,BLA203A30_S21_L002_R1_001.fastq.gz,,,
# WT_INPUT_REP3,BLA203A31_S21_L003_R1_001.fastq.gz,,,

read_csv(
    "senatorov-et-al-2023/GSE161948-PRJNA679976-SraRunTable.csv") %>%
    left_join(data.frame("file_name" = list.files(
        "/data/LGCP/freedman-chip/sra-files"),
        "full_file_name" = list.files(
            "/data/LGCP/freedman-chip/sra-files",
            full.names = T)) %>%
        mutate(Run = str_remove(file_name, ".fastq.gz") %>%
            str_remove("_[1|2]"))) %>%
    select(Run, file_name, ChIP_antibody, `Sample Name`, everything()) %>%
    as.data.frame() %>%
    left_join(acc_to_sample_df,
        by = c("GEO_Accession (exp)")) %>%
    select(sample, full_file_name) %>%
    mutate(pivot_names = if_else(str_detect(full_file_name, "_1.fastq.gz"),
        "fastq_1", "fastq_2")) %>%
    pivot_wider(names_from = pivot_names,
        values_from = full_file_name) %>%
    transmute(sample = sample,
        fastq_1 = fastq_1,
        fastq_2 = fastq_2,
        antibody = str_remove(sample, "LuCaP_.*_"),
        control = str_replace(sample, antibody, "INPUT")) %>%
    mutate(antibody = if_else(antibody == "INPUT",
        "",
        toupper(antibody)),
        control = if_else(antibody == "",
            "",
            control)) %>%
    arrange(desc(antibody)) %>%
    filter(!str_detect(sample, "_35_")) %>% # remove 35, no input available
    write_csv("senatorov-et-al-2023/design-GSE161948-PRJNA679976.csv")
