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

nextflow run nf-core/atacseq -r dev --input /data/capaldobj/lgcp/senatorov-et-al-2023/design-AR-GSE130408-PRJNA540151.csv \
-profile biowulf \
--aligner bwa \
--genome GRCh37 \
--igenomes_base 's3://ngi-igenomes/igenomes' \
--narrow_peak \
--read_length 75 \
--outdir '/data/LGCP/freedman-chip/lucap-only-ar-results/'
