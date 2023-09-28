acc_to_sample_df <- sample_sheet <- data.frame(
  stringsAsFactors = FALSE,
            Sample = c(1L,2L,3L,4L,5L,6L,7L,8L,
                       9L,10L,11L,12L,13L,14L,15L,16L,17L,18L,19L,
                       20L,21L,22L,23L,24L,25L,26L,27L,28L,29L,30L,31L,
                       32L,33L,34L),
                ID = c("M1","H1","170.2","M2","H2",
                       "170.3T","M1","H1","170.2","170.3T","M1","H1",
                       "170.2","M2","H2","170.3T","M1r2","H1r2","170.2r2",
                       "M2r2","H2r2","170.3r2","M1r2","H1r2","170.2r2",
                       "M2r2","H2r2","170.3r2","M1r2","H1r2","170.2r2","M2r2",
                       "H2r2","170.3r2"),
            Target = c("TE input","TE input",
                       "TE input","TE input","TE input","TE input","HNF1A","HNF1A",
                       "HNF1A","HNF1A","H3K27Ac","H3K27Ac","H3K27Ac",
                       "H3K27Ac","H3K27Ac","H3K27Ac","TE input","TE input",
                       "TE input","TE input","TE input","TE input","HNF1A",
                       "HNF1A","HNF1A","HNF1A","HNF1A","HNF1A","H3K27Ac",
                       "H3K27Ac","H3K27Ac","H3K27Ac","H3K27Ac","H3K27Ac"),
       Target.Type = c("input control",
                       "input control","input control","input control","input control",
                       "input control","TF negative control","TF",
                       "TF negative control","TF","H3K Acetylation","H3K Acetylation",
                       "H3K Acetylation","H3K Acetylation","H3K Acetylation",
                       "H3K Acetylation","input control","input control",
                       "input control","input control","input control",
                       "input control","TF negative control","TF","TF negative control",
                       "TF negative control","TF","TF","H3K Acetylation",
                       "H3K Acetylation","H3K Acetylation","H3K Acetylation",
                       "H3K Acetylation","H3K Acetylation"),
     Concentration = c(0.071777778,0.088,0.102444444,
                       0.114,0.098,0.222222222,0.0032,0.0603,0.0093,0.04,
                       0.125,0.124,0.256,0.19,0.217,0.283,39,40,82.5,
                       27,28,111.3,0.011,0.014,0.027,0.018,0.037,0.041,
                       0.27,0.208,0.14,0.226,0.51,0.4),
            Volume = c(30L,30L,30L,30L,30L,30L,
                       20L,20L,20L,20L,20L,20L,20L,20L,20L,20L,30L,30L,
                       30L,30L,30L,30L,20L,20L,20L,20L,20L,13L,20L,
                       20L,20L,20L,20L,20L)
)

library(tidyverse)
# GSE130408-PRJNA540151
acc_to_sample_df %>%
    mutate(sample_id = paste0("Sample_",
        Sample,
        "_Senatorov_",
        Sample)) %>%
    left_join(data.frame(
        "full_file_name" = list.files(
            "/data/capaldobj/incoming-nih-dme/CS035088-chip-seq-ilya-hnf1a/",
            full.names = T,
            recursive = T,
            pattern = ".fastq.gz$")) %>%
        mutate(sample_id = basename(dirname(full_file_name)),
            fastq_id = if_else(str_detect(full_file_name, "R2"),
                "fastq_2",
                "fastq_1"))) %>%
    pivot_wider(names_from = fastq_id,
        values_from = full_file_name) %>%
        unnest(cols = c(fastq_1, fastq_2)) %>%
    select(ID:Target.Type, sample_id:fastq_2) %>%
    transmute(sample = paste(ID, Target, sep = "_"),
        fastq_1 = fastq_1,
        fastq_2 = fastq_2,
        antibody = Target) %>%
    mutate(antibody = if_else(str_detect(antibody, "input"),
            "",
            antibody),
        control = str_replace(sample,
            "_HNF1A|_H3K27Ac",
            "_TE input"),
        control = if_else(antibody ==  "",
            "",
            control)) %>%
    as.data.frame() %>%
    arrange(desc(antibody)) %>%
    write_csv("senatorov-et-al-2023/design-hnf1a-chip.csv")
