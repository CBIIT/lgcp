library(tidyverse)

peak_files <- paste0("/data/capaldobj/incoming-nih-dme/",
    "CS035088-chip-seq-ilya-hnf1a/",
    "chipseq-grch37-results/bwa/mergedLibrary/macs2/narrowPeak/",
    "consensus/H3K27Ac/H3K27Ac.consensus_peaks.bed")

bam_files <- list.files(
    paste0("/data/capaldobj/incoming-nih-dme/CS035088-chip-seq-ilya-hnf1a/",
    "chipseq-grch37-results/bwa/mergedLibrary/"),
    pattern = "H3K27Ac.mLb.clN.sorted.bam$",
    full.names = T
)

for (sample_file in peak_files) {
    read_tsv(sample_file,
        col_names = F) %>%
    mutate(X1 = paste0("chr", X1)) %>%
    write_tsv(str_replace(sample_file,
        "consensus_peaks",
        "consensus_peaks.hg19"))
}

# rgt only supports hg19, so adding "chr" in front of all contigs is req
# run code below, then manually paste in samtools command
# samtools reheader -c "sed -e 's/SN:\([0-9XY]\+\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/'" h37.bam > hg19.bam
# swarm -f /data/capaldobj/lgcp/senatorov-et-al-2023/hnf1a-grch37-2-hg19-bam-converter.swarm -t 8 -g 64 --module samtools --partition ccr
paste0(bam_files,
    " > ",
    str_replace(bam_files, "mLb.clN.sorted.bam", "hg19.sorted.bam")) %>%
    data.frame(command = .) %>%
    write_csv("senatorov-et-al-2023/hnf1a-grch37-2-hg19-bam-converter.swarm",
        col_names = F)

# reindex bam files
# swarm -f lgcp/senatorov-et-al-2023/index-hg19-bam-converter.swarm -t 4 -g 16 --module samtools --partition ccr
paste0("samtools index ",
    str_replace(bam_files, "mLb.clN.sorted.bam", "hg19.sorted.bam")) %>%
    data.frame(command = .) %>%
    write_csv("senatorov-et-al-2023/index-hg19-bam-converter.swarm",
        col_names = F)

# swarm -f lgcp/senatorov-et-al-2023/rgt-footprinting.swarm -t 8 -g 64 --module rgt --partition ccr

paste0("rgt-hint footprinting ",
    "--organism=hg19 ",
    "--output-location=",
    "/data/capaldobj/incoming-nih-dme/CS035088-chip-seq-ilya-hnf1a/",
    "chipseq-grch37-results/Footprints ",
    "--output-prefix=",
    basename(bam_files) %>%
        str_remove(".mLb.clN.sorted.bam"),
    " --histone ",
    str_replace(bam_files, "mLb.clN.sorted.bam", "hg19.sorted.bam"),
    " ",
    str_replace(peak_files, "consensus_peaks", "consensus_peaks.hg19")) %>%
    data.frame(command = .) %>%
    write_csv("senatorov-et-al-2023/rgt-footprinting.swarm",
        col_names = F)

# swarm -f lgcp/senatorov-et-al-2023/rgt-motifmatching.swarm -t 8 -g 64 --module rgt --partition ccr

# rgt-motifanalysis matching --organism=hg19 --output-location=./MotifMatching --input-files ./Footprints/adeno_diff.bed

paste0("rgt-motifanalysis matching ",
    "--organism=hg19 ",
    "--output-location=",
    "/data/capaldobj/incoming-nih-dme/",
    "CS035088-chip-seq-ilya-hnf1a/chipseq-grch37-results/MotifMatching ",
    "--input-files ",
    list.files(
        paste0("/data/capaldobj/incoming-nih-dme/",
            "CS035088-chip-seq-ilya-hnf1a/chipseq-grch37-results/Footprints"),
        pattern = ".bed",
        full.names = T)) %>%
    data.frame(command = .) %>%
    write_csv("senatorov-et-al-2023/rgt-motifmatching.swarm",
        col_names = F)
