Incorporating Homer into LGCP ChIPSeq Pipeline
================
Brian Capaldo
3/14/2019

``` java
process createBigWig {
  tag "${bam.baseName - 'sortedByCoord.out'}"
  publishDir "${params.outdir}/bigwig", mode: 'copy'
  
  when:
    !params.skip_qc && !params.skip_genebody_coverage
  
  input:
    file bam from bam_for_genebody
  file index from bam_index_genebody
  
  output:
    file "*.bigwig" into bigwig_for_genebody
  
  script:
    """
  bamCoverage -b $bam -p ${task.cpus} -o ${bam.baseName}.bigwig
  """
}
```

``` bash
input_file=/data/khanlab/projects/ChIP_seq/DATA/Sample_RH4_H3K27ac_001_C_H5TLGBGXX/Sample_RH4_H3K27ac_001_C_H5TLGBGXX.dd.bam_7X_peaks.narrowPeak.bed
output_file=/data/khanlab/projects/ChIP_seq/DATA/Sample_RH4_H3K27ac_001_C_H5TLGBGXX/Sample_RH4_H3K27ac_001_C_H5TLGBGXX.dd.bam_7X_peaks.narrowPeak.bed.annotation.txt

annotatePeaks.pl $input_file hg19 > $output_file

#input_file=/data/khanlab/projects/ChIP_seq/data_by_file_type/test/DATA/Sample_RH4_H3K27ac_001_C_H5TLGBGXX/ROSE_out_narrow_1000/trimmed_AllEnhancers.sorted.table.txt_fpkm3.regular.bed
#output_dir=/data/khanlab/projects/ChIP_seq/data_by_file_type/test/DATA/Sample_RH4_H3K27ac_001_C_H5TLGBGXX/ROSE_out_narrow_1000/motif_regular
#preparse_dir=/data/khanlab/projects/ChIP_seq/data_by_file_type/test/scripts/preparse


findMotifsGenome.pl <peak/BED file> <genome> <output directory> -size 200 -mask -p <# threads> -preparse
  findMotifsGenome.pl $input_file $genome $output_dir -size $size -len $len -preparsedDir $preparsedDir -p $p

# Region Size ("-size <#>", "-size <#>,<#>", "-size given", default: 200)
# The size of the region used for motif finding is important.  If analyzing ChIP-Seq peaks from a transcription factor, Chuck would recommend 50 bp for establishing the primary motif bound by a given transcription factor and 200 bp for finding both primary and "co-enriched" motifs for a transcription factor.  When looking at histone marked regions, 500-1000 bp is probably a good idea (i.e. H3K4me or H3/H4 acetylated regions).  In theory, HOMER can work with very large regions (i.e. 10kb), but with the larger the regions comes more sequence and longer execution time.  These regions will be based off the center of the peaks.  If you prefer an offset, you can specify "-size -300,100" to search a region of size 400 that is centered 100 bp upstream of the peak center (useful if doing motif finding on putative TSS regions).  If you have variable length regions, use the option "-size given" and HOMER will use the exact regions that were used as input.
```
