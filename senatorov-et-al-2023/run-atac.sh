#! /bin/bash
#SBATCH --cpus-per-task=12
#SBATCH --partition=ccr
#SBATCH --mem=32gb
#SBATCH --time=3-00:00:00
#SBATCH --mail-type=BEGIN,TIME_LIMIT_90,END
#SBATCH --mail-user=capaldobj@nih.gov

module purge
module load nextflow
module load singularity
module load graphviz

nextflow run nf-core/atacseq --input design.csv \
-profile biowulf \
-resume \
--aligner bwa \
--genome GRCh37 \
--read_length 75 \
--narrow_peak \
--outdir '/data/LGCP/freedman-chip/results/'
