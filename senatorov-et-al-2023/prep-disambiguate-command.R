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
         sample_id, " -o ngs_disambiguated -a star ",
                          bam_file_human,
                          " ",
                          bam_file_mouse)) %>%
  select(command) %>%
  write_tsv("/data/capaldobj/lgcp/senatorov-et-al-2023/ngs_disambiguate.swarm", col_names = F)

# swarm -f /data/capaldobj/lgcp/senatorov-et-al-2023/ngs_disambiguate.swarm -t 8 -g 64 --partition ccr

# featureCounts \
#    -F SAF -O --fracOverlap 0.2 \
#    -p \
#    -T 6 \
#    -a H3K27Ac.consensus_peaks.saf \
#    -s 0 \
#    -o H3K27Ac.consensus_peaks.featureCounts.txt \
#    H2r2_H3K27Ac.mLb.clN.sorted.bam M1r2_H3K27Ac.mLb.clN.sorted.bam H2_H3K27Ac.mLb.clN.sorted.bam H1r2_H3K27Ac.mLb.clN.sorted.bam M2_H3K27Ac.mLb.clN.sorted.bam 170.3T_H3K27Ac.mLb.clN.sorted.bam 170.2r2_H3K27Ac.mLb.clN.sorted.bam H1_H3K27Ac.mLb.clN.sorted.bam 170.3r2_H3K27Ac.mLb.clN.sorted.bam M2r2_H3K27Ac.mLb.clN.sorted.bam 170.2_H3K27Ac.mLb.clN.sorted.bam M1_H3K27Ac.mLb.clN.sorted.bam

# /data/capaldobj/incoming-nih-dme/CS035088-chip-seq-ilya-hnf1a/chipseq-grch37-results/bwa/mergedLibrary/macs2/narrowPeak/consensus/HNF1A/HNF1A.consensus_peaks.saf

data.frame(bam_file = list.files(
  "/data/capaldobj/incoming-nih-dme/CS035088-chip-seq-ilya-hnf1a/ngs_disambiguated",
                                  pattern = ".disambiguatedSpeciesA.bam",
                                  full.names = T),
           sample_id = basename(
            list.files(
              "/data/capaldobj/incoming-nih-dme/CS035088-chip-seq-ilya-hnf1a/ngs_disambiguated/",
                                  pattern = ".disambiguatedSpeciesA.bam",
                                  full.names = T))) %>%
   mutate(sample_id = str_remove(sample_id, ".disambiguatedSpeciesA.bam"),
         out_file = paste0(
            "/data/capaldobj/incoming-nih-dme/CS035088-chip-seq-ilya-hnf1a/featureCounts/",
            sample_id,
            "_grch37_counts.txt"),
         command = paste0("featureCounts -T 4 -p -t exon -g gene_id -a /fdb/igenomes/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf -o ",
                          out_file,
                          " ",
                          bam_file)) %>%
  select(command) %>%
  write_tsv(
    "/data/capaldobj/lgcp/senatorov-et-al-2023/featureCounts-human.swarm",
    col_names = F)

data.frame(bam_file = list.files("/data/capaldobj/CS026803-rna-seq-yin-ivy-bone-marrow-metastatis/ngs_disambiguated/",
                                  pattern = ".disambiguatedSpeciesB.bam",
                                  full.names = T),
           sample_id = basename(list.files("/data/capaldobj/CS026803-rna-seq-yin-ivy-bone-marrow-metastatis/ngs_disambiguated/",
                                  pattern = ".disambiguatedSpeciesB.bam",
                                  full.names = T))) %>% 
   mutate(sample_id = str_remove(sample_id, ".disambiguatedSpeciesB.bam"),
         out_file = paste0("/data/capaldobj/CS026803-rna-seq-yin-ivy-bone-marrow-metastatis/featureCounts/", sample_id, "_grcm38_counts.txt"),
         command = paste0("featureCounts -T 4 -p -t exon -g gene_id -a /fdb/igenomes/Mus_musculus/Ensembl/GRCm38/Annotation/Genes/genes.gtf -o ",
                          out_file,
                          " ",
                          bam_file))%>% 
  select(command) %>% 
  write_tsv("/data/capaldobj/CS026803-rna-seq-yin-ivy-bone-marrow-metastatis/featureCounts-mouse.swarm", col_names = F)