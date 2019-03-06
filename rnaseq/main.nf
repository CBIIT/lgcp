#!/usr/bin/env nextflow
 
/*
 * Defines pipeline parameters in order to specify the refence genomes
 * and read pairs by using the command line options
 */
params.reads = "$baseDir/data/ggal/*_{1,2}.fq"
params.outdir = "/data/$USER/"
params.seqCenter = false
params.aligner = 'star'

clip_r1 = 3
clip_r2 = 0
three_prime_clip_r1 = 0
three_prime_clip_r2 = 3

wherearemyfiles = file("$baseDir/assets/where_are_my_files.txt")
biotypes_header = file("$baseDir/assets/biotypes_header.txt")


/*
 * Creates the `read_pairs` channel that emits for each read-pair a tuple containing
 * three elements: the pair ID, the first read-pair file and the second read-pair file
 */

Channel
        .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
        .into { raw_reads_fastqc; raw_reads_trimgalore }

process fastqc {
    tag "$name"
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set val(name), file(reads) from raw_reads_fastqc

    output:
    file "*_fastqc.{zip,html}" into fastqc_results

    script:
    """
    fastqc -q $reads
    """
}

/*
 * STEP 2 - Trim Galore!
 */
process trim_galore {
    tag "$name"
    publishDir "${params.outdir}/trim_galore", mode: 'copy',
	saveAs: {filename ->
            if (filename.indexOf("_fastqc") > 0) "FastQC/$filename"
            else if (filename.indexOf("trimming_report.txt") > 0) "logs/$filename"
            else if (!params.saveTrimmed && filename == "where_are_my_files.txt") filename
            else if (params.saveTrimmed && filename != "where_are_my_files.txt") filename
            else null
        }

    input:
    set val(name), file(reads) from raw_reads_trimgalore
    file wherearemyfiles

    output:
    file "*fq.gz" into trimmed_reads
    file "*trimming_report.txt" into trimgalore_results
    file "*_fastqc.{zip,html}" into trimgalore_fastqc_reports
    file "where_are_my_files.txt"


    script:
    c_r1 = clip_r1 > 0 ? "--clip_r1 ${clip_r1}" : ''
    c_r2 = clip_r2 > 0 ? "--clip_r2 ${clip_r2}" : ''
    tpc_r1 = three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${three_prime_clip_r1}" : ''
    tpc_r2 = three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${three_prime_clip_r2}" : ''
    """
    trim_galore --paired --fastqc --gzip $c_r1 $c_r2 $tpc_r1 $tpc_r2 $reads
    """
}
/*
 * STEP 3 - align with STAR
 */
// Function that checks the alignment rate of the STAR output
// and returns true if the alignment passed and otherwise false
skipped_poor_alignment = []
def check_log(logs) {
    def percent_aligned = 0;
    logs.eachLine { line ->
        if ((matcher = line =~ /Uniquely mapped reads %\s*\|\s*([\d\.]+)%/)) {
            percent_aligned = matcher[0][1]
        }
    }
    logname = logs.getBaseName() - 'Log.final'
    if(percent_aligned.toFloat() <= '5'.toFloat() ){
        log.info "#################### VERY POOR ALIGNMENT RATE! IGNORING FOR FURTHER DOWNSTREAM ANALYSIS! ($logname)    >> ${percent_aligned}% <<"
        skipped_poor_alignment << logname
        return false
    } else {
        log.info "          Passed alignment > star ($logname)   >> ${percent_aligned}% <<"
        return true
    }
}
if(params.aligner == 'star'){
    hisat_stdout = Channel.from(false)
    process star {
        tag "$prefix"
        publishDir "${params.outdir}/STAR", mode: 'copy',
            saveAs: {filename ->
                if (filename.indexOf(".bam") == -1) "logs/$filename"
                else if (!params.saveAlignedIntermediates && filename == "where_are_my_files.txt") filename
                else if (params.saveAlignedIntermediates && filename != "where_are_my_files.txt") filename
                else null
            }
        
        clusterOptions = '--partition ccr --cpus-per-task=16 --mem=128g --gres=lscratch:60'
        
        input:
        file reads from trimmed_reads
        file wherearemyfiles

        output:
        set file("*Log.final.out"), file ('*.bam') into star_aligned
        file "*.out" into alignment_logs
        file "*SJ.out.tab"
        file "*Log.out" into star_log
        file "where_are_my_files.txt"

        script:
        prefix = reads[0].toString() - ~/(_R1)?(_trimmed)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
        seqCenter = params.seqCenter ? "--outSAMattrRGline ID:$prefix 'CN:$params.seqCenter'" : ''
        """
        STAR --genomeDir "/data/capaldobj/genomes/starindex-ercc/" \
            --sjdbGTFfile "/data/capaldobj/genomes/grch37-ensembl-ercc.gtf" \
            --readFilesIn $reads  \
            --runThreadN 16 \
            --twopassMode Basic \
            --outWigType bedGraph \
            --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 60000000000 \
            --readFilesCommand zcat \
            --runDirPerm All_RWX \
            --outFileNamePrefix $prefix $seqCenter
        """
    }
    // Filter removes all 'aligned' channels that fail the check
    star_aligned
        .filter { logs, bams -> check_log(logs) }
        .flatMap {  logs, bams -> bams }
    .into { bam_count; bam_rseqc; bam_preseq; bam_markduplicates; bam_featurecounts; bam_stringtieFPKM; bam_for_genebody }
}

/*
 * STEP 8 Feature counts
 */
process featureCounts {
    tag "${bam_featurecounts.baseName - '.sorted'}"
    publishDir "${params.outdir}/featureCounts", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("biotype_counts") > 0) "biotype_counts/$filename"
            else if (filename.indexOf("_gene.featureCounts.txt.summary") > 0) "gene_count_summaries/$filename"
            else if (filename.indexOf("_gene.featureCounts.txt") > 0) "gene_counts/$filename"
            else "$filename"
        }

    clusterOptions = '--partition ccr --cpus-per-task=8 --mem=64g'
    
    input:
    file bam_featurecounts

    output:
    file "${bam_featurecounts.baseName}_gene.featureCounts.txt" into geneCounts, featureCounts_to_merge
    file "${bam_featurecounts.baseName}_gene.featureCounts.txt.summary" into featureCounts_logs

    script:
    // Try to get real sample name
    sample_name = bam_featurecounts.baseName - 'Aligned.sortedByCoord.out'
    """
    featureCounts -a "/data/capaldobj/genomes/grch37-ensembl-ercc.gtf" -g "gene_id" -F "GTF" -t "exon" -o ${bam_featurecounts.baseName}_gene.featureCounts.txt -p -C -B -T 8 $bam_featurecounts
    """
}
