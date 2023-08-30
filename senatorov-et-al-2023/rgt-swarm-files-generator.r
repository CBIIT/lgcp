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

for (sample_file in peak_files) {
    read_tsv(sample_file,
        col_names = F) %>%
    mutate(X1 = paste0("chr", X1)) %>%
    write_tsv(str_replace(sample_file, "peaks.narrowPeak", "peaks.hg19.narrowPeak"))
}

# rgt only supports hg19, so adding "chr" in front of all contigs is req
# run code below, then manually paste in samtools command
# samtools reheader -c "sed -e 's/SN:\([0-9XY]\+\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/'" h37.bam > hg19.bam
# swarm -f grch37-2-hg19-bam-converter.swarm -t 8 -g 64 --module samtools --partition ccr
paste0(bam_files,
    " > ",
    str_replace(bam_files, "sorted.bam", "hg19.sorted.bam")) %>%
    data.frame(command = .) %>%
    write_csv("senatorov-et-al-2023/grch37-2-hg19-bam-converter.swarm",
        col_names = F)

# swarm -f footprints.swarm -t 8 -g 64 --module rgt --partition ccr

paste0("rgt-hint footprinting ",
    "--organism=hg19 ",
    "--output-location=/data/LGCP/freedman-chip/lucap-only-k27ac-results/Footprints ",
    "--output-prefix=",
    basename(bam_files) %>%
        str_remove("_REP1.*$"),
    " --histone ",
    str_replace(bam_files, "sorted.bam", "hg19.sorted.bam"),
    " ",
    str_replace(peak_files, "peaks.narrowPeak", "peaks.hg19.narrowPeak")) %>%
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