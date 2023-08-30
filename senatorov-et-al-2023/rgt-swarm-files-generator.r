library(tidyverse)

peak_files <- list.files(
    "/data/LGCP/freedman-chip/lucap-only-k27ac-results/bwa/merged_library/macs2/narrow_peak",
    pattern = "peaks.narrowPeak",
    full.names = T)

bam_files <- list.files(
    "/data/LGCP/freedman-chip/lucap-only-k27ac-results/bwa/merged_library",
    pattern = "sorted.bam$",
    full.names = T
)

# swarm -f footprints.swarm -t 8 -g 64 --module rgt --partition ccr

paste0("rgt-hint footprinting ",
    "--organism=hg19 ",
    "--output-location=/data/LGCP/freedman-chip/lucap-only-k27ac-results/Footprints ",
    "--output-prefix=",
    basename(bam_files) %>%
        str_remove("_REP1.*$"),
    " --histone ",
    bam_files,
    " ",
    peak_files) %>%
    data.frame(command = .) %>%
    write_csv("senatorov-et-al-2023/rgt-footprinting.swarm",
        col_names = F)

# swarm -f motif-matching.swarm -t 8 -g 64 --module rgt --partition ccr

# rgt-motifanalysis matching --organism=hg19 --output-location=./MotifMatching --input-files ./Footprints/adeno_diff.bed

paste0("rgt-motifanalysis matching ",
    "--organism=hg19 ",
    "--output-location=/data/LGCP/freedman-chip/lucap-only-k27ac-results/MotifMatching ",
    "--input-files ",
    list.files("/data/LGCP/freedman-chip/lucap-only-k27ac-results/Footprints/")) %>%
    data.frame(command = .) %>%
    write_csv("senatorov-et-al-2023/rgt-motifmatching.swarm",
        col_names = F)