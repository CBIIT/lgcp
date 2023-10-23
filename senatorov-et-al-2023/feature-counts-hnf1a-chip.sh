#!/bin/bash
#SBATCH -J featurecountshijack
#SBATCH -c 6
#SBATCH -t 08:00:00
#SBATCH --mem 36864M
#SBATCH  --gres=lscratch:200 --partition ccr 

module load subread

featureCounts \
    -F SAF -O --fracOverlap 0.2 \
    -p \
    -T 6 \
    -a /data/capaldobj/incoming-nih-dme/CS035088-chip-seq-ilya-hnf1a/chipseq-grch37-results/bwa/mergedLibrary/macs2/narrowPeak/consensus/H3K27Ac/H3K27Ac.consensus_peaks.saf \
    -s 0 \
    -o H3K27Ac.grch37.consensus_peaks.featureCounts.txt \
    H2r2_H3K27Ac.disambiguatedSpeciesA.bam M1r2_H3K27Ac.disambiguatedSpeciesA.bam H2_H3K27Ac.disambiguatedSpeciesA.bam H1r2_H3K27Ac.disambiguatedSpeciesA.bam M2_H3K27Ac.disambiguatedSpeciesA.bam 170.3T_H3K27Ac.disambiguatedSpeciesA.bam 170.2r2_H3K27Ac.disambiguatedSpeciesA.bam H1_H3K27Ac.disambiguatedSpeciesA.bam 170.3r2_H3K27Ac.disambiguatedSpeciesA.bam M2r2_H3K27Ac.disambiguatedSpeciesA.bam 170.2_H3K27Ac.disambiguatedSpeciesA.bam M1_H3K27Ac.disambiguatedSpeciesA.bam

featureCounts \
    -F SAF -O --fracOverlap 0.2 \
    -p \
    -T 6 \
    -a /data/capaldobj/incoming-nih-dme/CS035088-chip-seq-ilya-hnf1a/chipseq-grch37-results/bwa/mergedLibrary/macs2/narrowPeak/consensus/HNF1A/HNF1A.consensus_peaks.saf \
    -s 0 \
    -o HNF1A.grch37.consensus_peaks.featureCounts.txt \
    170.2_HNF1A.disambiguatedSpeciesA.bam 170.2r2_HNF1A.disambiguatedSpeciesA.bam 170.3r2_HNF1A.disambiguatedSpeciesA.bam 170.3T_HNF1A.disambiguatedSpeciesA.bam H1_HNF1A.disambiguatedSpeciesA.bam H1r2_HNF1A.disambiguatedSpeciesA.bam H2r2_HNF1A.disambiguatedSpeciesA.bam M1_HNF1A.disambiguatedSpeciesA.bam M1r2_HNF1A.disambiguatedSpeciesA.bam M2r2_HNF1A.disambiguatedSpeciesA.bam

featureCounts \
    -F SAF -O --fracOverlap 0.2 \
    -p \
    -T 6 \
    -a /data/capaldobj/incoming-nih-dme/CS035088-chip-seq-ilya-hnf1a/chipseq-results/bwa/mergedLibrary/macs2/narrowPeak/consensus/H3K27Ac/H3K27Ac.consensus_peaks.saf \
    -s 0 \
    -o H3K27Ac.grcm38.consensus_peaks.featureCounts.txt \
    170.2_H3K27Ac.disambiguatedSpeciesB.bam 170.2r2_H3K27Ac.disambiguatedSpeciesB.bam 170.3r2_H3K27Ac.disambiguatedSpeciesB.bam 170.3T_H3K27Ac.disambiguatedSpeciesB.bam H1_H3K27Ac.disambiguatedSpeciesB.bam H1r2_H3K27Ac.disambiguatedSpeciesB.bam H2_H3K27Ac.disambiguatedSpeciesB.bam H2r2_H3K27Ac.disambiguatedSpeciesB.bam M1_H3K27Ac.disambiguatedSpeciesB.bam M1r2_H3K27Ac.disambiguatedSpeciesB.bam M2_H3K27Ac.disambiguatedSpeciesB.bam M2r2_H3K27Ac.disambiguatedSpeciesB.bam

featureCounts \
    -F SAF -O --fracOverlap 0.2 \
    -p \
    -T 6 \
    -a /data/capaldobj/incoming-nih-dme/CS035088-chip-seq-ilya-hnf1a/chipseq-results/bwa/mergedLibrary/macs2/narrowPeak/consensus/HNF1A/HNF1A.consensus_peaks.saf \
    -s 0 \
    -o HNF1A.grcm38.consensus_peaks.featureCounts.txt \
    170.2_HNF1A.disambiguatedSpeciesB.bam 170.2r2_HNF1A.disambiguatedSpeciesB.bam 170.3r2_HNF1A.disambiguatedSpeciesB.bam 170.3T_HNF1A.disambiguatedSpeciesB.bam H1_HNF1A.disambiguatedSpeciesB.bam H1r2_HNF1A.disambiguatedSpeciesB.bam H2r2_HNF1A.disambiguatedSpeciesB.bam M1_HNF1A.disambiguatedSpeciesB.bam M1r2_HNF1A.disambiguatedSpeciesB.bam M2r2_HNF1A.disambiguatedSpeciesB.bam
