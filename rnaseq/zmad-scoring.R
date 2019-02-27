library(tidyverse)
library(edgeR)
library(annotables)
library(readxl)
library(clusterProfiler)

count_data <- read_csv("")

# create tidy data frame of log expression values
annotated_cpm_df <- count_data %>%
  gather(sample_id,
         cpm_tmm,
         -c(ensgene:description)) %>%
  mutate(cpm_tmm = as.numeric(cpm_tmm),
    cpm_tmm_raw = 2^(cpm_tmm)) %>%
  filter(cpm_tmm_raw > 2,
         !is.na(cpm_tmm))

# calculate mean absolute deviation for each gene across all samples
median_absolute_deviations <- annotated_cpm_df %>%
  group_by(ensgene) %>%
  summarise(median_abs_dev = mad(cpm_tmm),
            median = median(cpm_tmm))
# calculate modified Z-score by multipling the difference of the cpm_tmm and mean for each gene times 0.6745 (just because), and dividing by the mean absolute deviation
modified_z_score <- annotated_cpm_df %>%
  left_join(median_absolute_deviations) %>%
  mutate(mod_z = (0.6745*(cpm_tmm - median))/median_abs_dev)

modified_z_score <- modified_z_score %>%
  separate(ensgene,
           c("ensgene", "version")) %>%
  left_join(grch37 %>%
              dplyr::select(entrez, symbol, ensgene))

# create function to calculate scores
really_big_df <- updated_signatures_df %>%
  left_join(modified_z_score,
            by = c("gene" = "symbol")) %>%
  group_by(ont, sample_id) %>%
  summarise(ont_score = mean(mod_z)) 