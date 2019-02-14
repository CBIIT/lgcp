#!/usr/bin/env nextflow

/*
========================================================================================
                  C H I P - S E Q   B E S T   P R A C T I C E
========================================================================================
 ChIP-seq Best Practice Analysis Pipeline. Started May 2016.
 #### Homepage / Documentation
 https://github.com/nf-core/chipseq
 @#### Authors
 Chuan Wang <chuan.wang@scilifelab.se>
 Phil Ewels <phil.ewels@scilifelab.se>
 Alex Peltzer <alexander.peltzer@qbic.uni-tuebingen.de>
----------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------
Pipeline overview:
 - 1:   FastQC for raw sequencing reads quality control
 - 2:   Trim Galore! for adapter trimming
 - 3.1: BWA alignment against reference genome
 - 3.2: Post-alignment processing and format conversion
 - 3.3: Statistics about mapped reads
 - 4:   Picard for duplicate read identification
 - 5:   Statistics about read counts
 - 6.1: Phantompeakqualtools for normalized strand cross-correlation (NSC) and relative strand cross-correlation (RSC)
 - 6.2: Summarize NSC and RSC
 - 7:   deepTools for fingerprint, coverage bigwig, and correlation plots of reads over genome-wide bins
 - 8:   NGSplot for distribution of reads around transcription start sites (TSS) and gene bodies
 - 9.1: MACS for peak calling
 - 9.2: Saturation analysis using MACS when specified
 - 10:  Post peak calling processing: blacklist filtering and annotation
 - 11:  MultiQC
 - 12:  Output Description HTML
 ----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info"""
    =========================================
     nf-core/chipseq : ChIP-Seq Best Practice v${params.version}
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/chipseq --reads '*_R{1,2}.fastq.gz' --genome GRCh37 --macsconfig 'macssetup.config' -profile uppmax

    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes).
      --genome                      Name of iGenomes reference
      --macsconfig                  Configuration file for peaking calling using MACS. Format: ChIPSampleID,CtrlSampleID,AnalysisID
      -profile                      Hardware config to use. uppmax / uppmax_modules / docker / aws

    Options:
      --singleEnd                   Specifies that the input is single end reads
      --allow_multi_align           Secondary alignments and unmapped reads are also reported in addition to primary alignments
      --saturation                  Run saturation analysis by peak calling with subsets of reads
      --broad                       Run MACS with the --broad flag
      --blacklist_filtering         Filter ENCODE blacklisted regions from ChIP-seq peaks. It only works when --genome is set as GRCh37 or GRCm38

    Presets:
      --extendReadsLen [int]        Number of base pairs to extend the reads for the deepTools analysis. Default: 100

    References
      --fasta                       Path to Fasta reference
      --bwa_index                   Path to BWA index
      --gtf                         Path to GTF file (Ensembl format)
      --blacklist                   Path to blacklist regions (.BED format), used for filtering out called peaks. Note that --blacklist_filtering is required
      --saveReference               Save the generated reference files in the Results directory.
      --saveAlignedIntermediates    Save the intermediate BAM files from the Alignment step  - not done by default

    Trimming options
      --notrim                      Specifying --notrim will skip the adapter trimming step.
      --saveTrimmed                 Save the trimmed Fastq files in the the Results directory.
      --clip_r1 [int]               Instructs Trim Galore to remove bp from the 5' end of read 1 (or single-end reads)
      --clip_r2 [int]               Instructs Trim Galore to remove bp from the 5' end of read 2 (paired-end reads only)
      --three_prime_clip_r1 [int]   Instructs Trim Galore to remove bp from the 3' end of read 1 AFTER adapter/quality trimming has been performed
      --three_prime_clip_r2 [int]   Instructs Trim Galore to re move bp from the 3' end of read 2 AFTER adapter/quality trimming has been performed

    Other options:
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      --rlocation                   Location to save R-libraries used in the pipeline. Default value is ~/R/nxtflow_libs/
      --clusterOptions              Extra SLURM options, used in conjunction with Uppmax.config
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help emssage
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

// Configurable variables
params.name = false
params.project = false
params.genome = false
params.genomes = false
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.bwa_index = params.genome ? params.genomes[ params.genome ].bwa ?: false : false
params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
params.blacklist = params.genome ? params.genomes[ params.genome ].blacklist ?: false : false

// R library locations
params.rlocation = false
if (params.rlocation){
    nxtflow_libs = file(params.rlocation)
    nxtflow_libs.mkdirs()
}

multiqc_config = file(params.multiqc_config)
output_docs = file("$baseDir/docs/output.md")

// Custom trimming options
params.clip_r1 = 0
params.clip_r2 = 0
params.three_prime_clip_r1 = 0
params.three_prime_clip_r2 = 0

// Validate inputs
macsconfig = file(params.macsconfig)
if( !macsconfig.exists() ) exit 1, "Missing MACS config: '$macsconfig'. Specify path with --macsconfig"
if( params.bwa_index ){
    bwa_index = Channel
        .fromPath(params.bwa_index)
        .ifEmpty { exit 1, "BWA index not found: ${params.bwa_index}" }
} else if ( params.fasta ){
    fasta = file(params.fasta)
    if( !fasta.exists() ) exit 1, "Fasta file not found: ${params.fasta}"
} else {
    exit 1, "No reference genome specified!"
}
gtf = false
if( params.gtf ){
    gtf = file(params.gtf)
    if( !gtf.exists() ) exit 1, "GTF file not found: ${params.gtf}."
}
if ( params.blacklist_filtering ){
    blacklist = file(params.blacklist)
    if( !blacklist.exists() ) exit 1, "Blacklist file not found: ${params.blacklist}"
}
if( workflow.profile == 'standard' && !params.project ) exit 1, "No UPPMAX project ID found! Use --project"

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

if(params.readPaths){
    if(params.singleEnd){
        Channel
            .from(params.readPaths)
            .map { row -> [ row[0], [file(row[1][0])]] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            .into { raw_reads_fastqc; raw_reads_trimgalore }
    } else {
        Channel
            .from(params.readPaths)
            .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            .into { raw_reads_fastqc; raw_reads_trimgalore }
    }
} else {
    Channel
        .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }
        .into { raw_reads_fastqc; raw_reads_trimgalore }
}


/*
 * Create a channel for macs config file
 */
Channel
    .from(macsconfig.readLines())
    .map { line ->
        list = line.split(',')
        chip_sample_id = list[0]
        ctrl_sample_id = list[1]
        analysis_id = list[2]
        [ chip_sample_id, ctrl_sample_id, analysis_id ]
    }
    .into{ macs_para; saturation_para }

/*
 * Reference to use for MACS, ngs.plot.r and annotation
 */
def REF_macs = false
def REF_ngsplot = false
if (params.genome == 'GRCh37'){ REF_macs = 'hs'; REF_ngsplot = 'hg19' }
else if (params.genome == 'GRCm38'){ REF_macs = 'mm'; REF_ngsplot = 'mm10' }


log.info """=======================================================
                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\'
    |\\ | |__  __ /  ` /  \\ |__) |__         }  {
    | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                          `._,._,\'

 nf-core/chipseq : ChIP-Seq Best Practice v${params.version}
======================================================="""
def summary = [:]
summary['Run Name']            = custom_runName ?: workflow.runName
summary['Reads']               = params.reads
summary['Data Type']           = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Genome']              = params.genome
if(params.bwa_index)  summary['BWA Index'] = params.bwa_index
else if(params.fasta) summary['Fasta Ref'] = params.fasta
if(params.gtf)  summary['GTF File'] = params.gtf
summary['Multiple alignments'] = params.allow_multi_align
summary['MACS Config']         = params.macsconfig
summary['Saturation analysis'] = params.saturation
summary['MACS broad peaks']    = params.broad
summary['Blacklist filtering'] = params.blacklist_filtering
if( params.blacklist_filtering ) summary['Blacklist BED'] = params.blacklist
summary['Extend Reads']        = "$params.extendReadsLen bp"
summary['Container']           = workflow.container
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Current home']        = "$HOME"
summary['Current user']        = "$USER"
summary['Current path']        = "$PWD"
summary['Working dir']         = workflow.workDir
summary['Output dir']          = params.outdir
summary['R libraries']         = params.rlocation
summary['Script dir']          = workflow.projectDir
summary['Save Reference']      = params.saveReference
summary['Save Trimmed']        = params.saveTrimmed
summary['Save Intermeds']      = params.saveAlignedIntermediates
if( params.notrim ){
    summary['Trimming Step'] = 'Skipped'
} else {
    summary['Trim R1'] = params.clip_r1
    summary['Trim R2'] = params.clip_r2
    summary["Trim 3' R1"] = params.three_prime_clip_r1
    summary["Trim 3' R2"] = params.three_prime_clip_r2
}
summary['Config Profile'] = workflow.profile
if(params.project) summary['UPPMAX Project'] = params.project
if(params.email) summary['E-mail Address'] = params.email
if(workflow.commitId) summary['Pipeline Commit']= workflow.commitId
log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "===================================="

// Check that Nextflow version is up to date enough
// try / throw / catch works for NF versions < 0.30.1 when this was implemented
nf_required_version = '0.30.1'
try {
    if( ! nextflow.version.matches(">= $nf_required_version") ){
        throw GroovyException('Nextflow version too old')
    }
} catch (all) {
    log.error "====================================================\n" +
              "  Nextflow version $nf_required_version required! You are running v$workflow.nextflow.version.\n" +
              "  Pipeline execution will continue, but things may break.\n" +
              "  Please run `nextflow self-update` to update Nextflow.\n" +
              "============================================================"
}

// Show a big error message if we're running on the base config and an uppmax cluster
if( workflow.profile == 'standard'){
    if ( "hostname".execute().text.contains('.uppmax.uu.se') ) {
        log.error "====================================================\n" +
                  "  WARNING! You are running with the default 'standard'\n" +
                  "  pipeline config profile, which runs on the head node\n" +
                  "  and assumes all software is on the PATH.\n" +
                  "  ALL JOBS ARE RUNNING LOCALLY and stuff will probably break.\n" +
                  "  Please use `-profile uppmax` to run on UPPMAX clusters.\n" +
                  "============================================================"
    }
}

// Show a big warning message if we're not running MACS
if (!REF_macs){
    def warnstring = params.genome ? "Reference '${params.genome}' not supported by" : 'No reference supplied for'
    log.warn "=======================================================\n" +
             "  WARNING! $warnstring MACS, ngs_plot\n" +
             "  and annotation. Steps for MACS, ngs_plot and annotation\n" +
             "  will be skipped. Use '--genome GRCh37' or '--genome GRCm38'\n" +
             "  to run these steps.\n" +
             "==============================================================="
}


/*
 * PREPROCESSING - Build BWA index
 */
if(!params.bwa_index && fasta){
    process makeBWAindex {
        tag fasta
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file fasta from fasta

        output:
        file "BWAIndex" into bwa_index

        script:
        """
        bwa index -a bwtsw $fasta
        mkdir BWAIndex && mv ${fasta}* BWAIndex
        """
    }
}


/*
 * STEP 1 - FastQC
 */
process fastqc {
    tag "$name"
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set val(name), file(reads) from raw_reads_fastqc

    output:
    file '*_fastqc.{zip,html}' into fastqc_results
    file '.command.out' into fastqc_stdout

    script:
    """
    fastqc -q $reads
    """
}


/*
 * STEP 2 - Trim Galore!
 */
if(params.notrim){
    trimmed_reads = read_files_trimming
    trimgalore_results = []
    trimgalore_fastqc_reports = []
} else {
    process trim_galore {
        tag "$name"
        publishDir "${params.outdir}/trim_galore", mode: 'copy',
            saveAs: {filename ->
                if (filename.indexOf("_fastqc") > 0) "FastQC/$filename"
                else if (filename.indexOf("trimming_report.txt") > 0) "logs/$filename"
                else params.saveTrimmed ? filename : null
            }

        input:
        set val(name), file(reads) from raw_reads_trimgalore

        output:
        file '*.fq.gz' into trimmed_reads
        file '*trimming_report.txt' into trimgalore_results
        file "*_fastqc.{zip,html}" into trimgalore_fastqc_reports

        script:
        c_r1 = params.clip_r1 > 0 ? "--clip_r1 ${params.clip_r1}" : ''
        c_r2 = params.clip_r2 > 0 ? "--clip_r2 ${params.clip_r2}" : ''
        tpc_r1 = params.three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${params.three_prime_clip_r1}" : ''
        tpc_r2 = params.three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${params.three_prime_clip_r2}" : ''
        if (params.singleEnd) {
            """
            trim_galore --fastqc --gzip $c_r1 $tpc_r1 $reads
            """
        } else {
            """
            trim_galore --paired --fastqc --gzip $c_r1 $c_r2 $tpc_r1 $tpc_r2 $reads
            """
        }
    }
}


/*
 * STEP 3.1 - align with bwa
 */
process bwa {
    tag "$prefix"
    publishDir path: { params.saveAlignedIntermediates ? "${params.outdir}/bwa" : params.outdir }, mode: 'copy',
               saveAs: {filename -> params.saveAlignedIntermediates ? filename : null }

    input:
    file reads from trimmed_reads
    file index from bwa_index.first()

    output:
    file '*.bam' into bwa_bam

    script:
    prefix = reads[0].toString() - ~/(.R1)?(_1)?(_R1)?(_trimmed)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
    filtering = params.allow_multi_align ? '' : "| samtools view -b -q 1 -F 4 -F 256"
    """
    bwa mem -M ${index}/genome.fa $reads | samtools view -bT $index - $filtering > ${prefix}.bam
    """
}


/*
 * STEP 3.2 - post-alignment processing
 */

process samtools {
    tag "${bam.baseName}"
    publishDir path: "${params.outdir}/bwa", mode: 'copy',
               saveAs: { filename ->
                   if (filename.indexOf(".stats.txt") > 0) "stats/$filename"
                   else params.saveAlignedIntermediates ? filename : null
               }

    input:
    file bam from bwa_bam

    output:
    file '*.sorted.bam' into bam_picard, bam_for_mapped
    file '*.sorted.bam.bai' into bwa_bai, bai_for_mapped
    file '*.sorted.bed' into bed_total
    file '*.stats.txt' into samtools_stats

    script:
    """
    samtools sort $bam -o ${bam.baseName}.sorted.bam
    samtools index ${bam.baseName}.sorted.bam
    bedtools bamtobed -i ${bam.baseName}.sorted.bam | sort -k 1,1 -k 2,2n -k 3,3n -k 6,6 > ${bam.baseName}.sorted.bed
    samtools stats ${bam.baseName}.sorted.bam > ${bam.baseName}.stats.txt
    """
}


/*
 * STEP 3.3 - Statistics about mapped and unmapped reads against ref genome
 */

process bwa_mapped {
    tag "${input_files[0].baseName}"
    publishDir "${params.outdir}/bwa/mapped", mode: 'copy'

    input:
    file input_files from bam_for_mapped.collect()
    file bai from bai_for_mapped.collect()

    output:
    file 'mapped_refgenome.txt' into bwa_mapped

    script:
    """
    for i in $input_files
    do
      samtools idxstats \${i} | awk -v filename="\${i}" '{mapped+=\$3; unmapped+=\$4} END {print filename,"\t",mapped,"\t",unmapped}'
    done > mapped_refgenome.txt
    """
}

/*
 * STEP 4 Picard
 */

process picard {
    tag "$prefix"
    publishDir "${params.outdir}/picard", mode: 'copy'

    input:
    file bam from bam_picard

    output:
    file '*.dedup.sorted.bam' into bam_dedup_spp, bam_dedup_ngsplot, bam_dedup_deepTools, bam_dedup_macs, bam_dedup_saturation
    file '*.dedup.sorted.bam.bai' into bai_dedup_deepTools, bai_dedup_ngsplot, bai_dedup_macs, bai_dedup_saturation
    file '*.dedup.sorted.bed' into bed_dedup
    file '*.picardDupMetrics.txt' into picard_reports

    script:
    prefix = bam[0].toString() - ~/(\.sorted)?(\.bam)?$/
    if( task.memory == null ){
        log.warn "[Picard MarkDuplicates] Available memory not known - defaulting to 6GB ($prefix)"
        avail_mem = 6000
    } else {
        avail_mem = task.memory.toMega()
        if( avail_mem <= 0){
            avail_mem = 6000
            log.warn "[Picard MarkDuplicates] Available memory 0 - defaulting to 6GB ($prefix)"
        } else if( avail_mem < 250){
            avail_mem = 250
            log.warn "[Picard MarkDuplicates] Available memory under 250MB - defaulting to 250MB ($prefix)"
        }
    }
    """
    picard MarkDuplicates \\
        INPUT=$bam \\
        OUTPUT=${prefix}.dedup.bam \\
        ASSUME_SORTED=true \\
        REMOVE_DUPLICATES=true \\
        METRICS_FILE=${prefix}.picardDupMetrics.txt \\
        VALIDATION_STRINGENCY=LENIENT \\
        PROGRAM_RECORD_ID='null'

    samtools sort ${prefix}.dedup.bam -o ${prefix}.dedup.sorted.bam
    samtools index ${prefix}.dedup.sorted.bam
    bedtools bamtobed -i ${prefix}.dedup.sorted.bam | sort -k 1,1 -k 2,2n -k 3,3n -k 6,6 > ${prefix}.dedup.sorted.bed
    """
}


/*
 * STEP 5 Read_count_statistics
 */

process countstat {
    tag "${input[0].baseName}"
    publishDir "${params.outdir}/countstat", mode: 'copy'

    input:
    file input from bed_total.mix(bed_dedup).toSortedList()

    output:
    file 'read_count_statistics.txt' into countstat_results

    script:
    """
    countstat.pl $input
    """
}


/*
 * STEP 6.1 Phantompeakqualtools
 * TODO: The "run_spp.R" script is still missing here!
 */

process phantompeakqualtools {
    tag "$prefix"
    publishDir "${params.outdir}/phantompeakqualtools", mode: 'copy',
                saveAs: {filename -> filename.indexOf(".out") > 0 ? "logs/$filename" : "$filename"}

    input:
    file bam from bam_dedup_spp

    output:
    file '*.pdf' into spp_results
    file '*.spp.out' into spp_out, spp_out_mqc

    script:
    prefix = bam[0].toString() - ~/(\.dedup)?(\.sorted)?(\.bam)?$/
    """
    run_spp.r -c="$bam" -savp -out="${prefix}.spp.out"
    """
}


/*
 * STEP 6.2 Combine and calculate NSC & RSC
 */

process calculateNSCRSC {
    tag "${spp_out_list[0].baseName}"
    publishDir "${params.outdir}/phantompeakqualtools", mode: 'copy'

    input:
    file spp_out_list from spp_out.collect()

    output:
    file 'cross_correlation_processed.txt' into calculateNSCRSC_results

    script:
    """
    cat $spp_out_list > cross_correlation.txt
    calculateNSCRSC.r cross_correlation.txt
    """
}


/*
 * STEP 7 deepTools
 */

process deepTools {
    tag "${bam[0].baseName}"
    publishDir "${params.outdir}/deepTools", mode: 'copy'

    input:
    file bam from bam_dedup_deepTools.collect()
    file bai from bai_dedup_deepTools.collect()

    output:
    file '*.{txt,pdf,png,npz,bw}' into deepTools_results
    file '*.txt' into deepTools_multiqc

    script:
    if (!params.singleEnd) {
        """
        bamPEFragmentSize \\
            --binSize 1000 \\
            --bamfiles $bam \\
            --histogram fragment_length_distribution_histogram.png \\
            --plotTitle "Fragment Length Distribution"
        """
    }
    if(bam instanceof Path){
        log.warn("Only 1 BAM file - skipping multiBam deepTool steps")
        """
        plotFingerprint \\
            -b $bam \\
            --plotFile ${bam.baseName}_fingerprints.pdf \\
            --outRawCounts ${bam.baseName}_fingerprint.txt \\
            --extendReads ${params.extendReadsLen} \\
            --skipZeros \\
            --ignoreDuplicates \\
            --numberOfSamples 50000 \\
            --binSize 500 \\
            --plotFileFormat pdf \\
            --plotTitle "${bam.baseName} Fingerprints"

        bamCoverage \\
           -b $bam \\
           --extendReads ${params.extendReadsLen} \\
           --normalizeUsing RPKM \\
           -o ${bam}.bw
        """
    } else {
        """
        plotFingerprint \\
            -b $bam \\
            --plotFile fingerprints.pdf \\
            --outRawCounts fingerprint.txt \\
            --extendReads ${params.extendReadsLen} \\
            --skipZeros \\
            --ignoreDuplicates \\
            --numberOfSamples 50000 \\
            --binSize 500 \\
            --plotFileFormat pdf \\
            --plotTitle "Fingerprints"

        for bamfile in ${bam}
        do
            bamCoverage \\
              -b \$bamfile \\
              --extendReads ${params.extendReadsLen} \\
              --normalizeUsing RPKM \\
              -o \${bamfile}.bw
        done

        multiBamSummary \\
            bins \\
            --binSize 10000 \\
            --bamfiles $bam \\
            -out multiBamSummary.npz \\
            --extendReads ${params.extendReadsLen} \\
            --ignoreDuplicates \\
            --centerReads

        plotCorrelation \\
            -in multiBamSummary.npz \\
            -o scatterplot_PearsonCorr_multiBamSummary.png \\
            --outFileCorMatrix scatterplot_PearsonCorr_multiBamSummary.txt \\
            --corMethod pearson \\
            --skipZeros \\
            --removeOutliers \\
            --plotTitle "Pearson Correlation of Read Counts" \\
            --whatToPlot scatterplot

        plotCorrelation \\
            -in multiBamSummary.npz \\
            -o heatmap_SpearmanCorr_multiBamSummary.png \\
            --outFileCorMatrix heatmap_SpearmanCorr_multiBamSummary.txt \\
            --corMethod spearman \\
            --skipZeros \\
            --plotTitle "Spearman Correlation of Read Counts" \\
            --whatToPlot heatmap \\
            --colorMap RdYlBu \\
            --plotNumbers

        plotPCA \\
            -in multiBamSummary.npz \\
            -o pcaplot_multiBamSummary.png \\
            --plotTitle "Principal Component Analysis Plot" \\
            --outFileNameData pcaplot_multiBamSummary.txt
        """
    }
}


/*
 * STEP 8 Ngsplot
 * TODO ngs.plot.R is missing too!
 */

process ngsplot {
    tag "${input_bam_files[0].baseName}"
    publishDir "${params.outdir}/ngsplot", mode: 'copy'

    input:
    file input_bam_files from bam_dedup_ngsplot.collect()
    file input_bai_files from bai_dedup_ngsplot.collect()

    output:
    file '*.pdf' into ngsplot_results

    when: REF_ngsplot

    script:
    """
    ngs_config_generate.r $input_bam_files

    ngs.plot.r \\
        -G $REF_ngsplot \\
        -R genebody \\
        -C ngsplot_config \\
        -O Genebody \\
        -D ensembl \\
        -FL 300

    ngs.plot.r \\
        -G $REF_ngsplot \\
        -R tss \\
        -C ngsplot_config \\
        -O TSS \\
        -FL 300
    """
}


/*
 * STEP 9.1 MACS
 */

process macs {
    tag "${bam_for_macs[0].baseName}"
    publishDir "${params.outdir}/macs", mode: 'copy'

    input:
    file bam_for_macs from bam_dedup_macs.collect()
    file bai_for_macs from bai_dedup_macs.collect()
    set chip_sample_id, ctrl_sample_id, analysis_id from macs_para

    output:
    file '*.{bed,r,narrowPeak}' into macs_results
    file '*.xls' into macs_peaks

    when: REF_macs

    script:
    def ctrl = ctrl_sample_id == '' ? '' : "-c ${ctrl_sample_id}.dedup.sorted.bam"
    broad = params.broad ? "--broad" : ''
    """
    macs2 callpeak \\
        -t ${chip_sample_id}.dedup.sorted.bam \\
        $ctrl \\
        $broad \\
        -f BAM \\
        -g $REF_macs \\
        -n $analysis_id \\
        -q 0.01
    """
}


/*
 * STEP 9.2 Saturation analysis
 */
if (params.saturation) {

  process saturation {
     tag "${bam_for_saturation[0].baseName}"
     publishDir "${params.outdir}/macs/saturation", mode: 'copy'

     input:
     file bam_for_saturation from bam_dedup_saturation.collect()
     file bai_for_saturation from bai_dedup_saturation.collect()
     set chip_sample_id, ctrl_sample_id, analysis_id from saturation_para
     each sampling from 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0

     output:
     file '*.xls' into saturation_results

     when: REF_macs

     script:
     def ctrl = ctrl_sample_id == '' ? '' : "-c ${ctrl_sample_id}.dedup.sorted.bam"
     broad = params.broad ? "--broad" : ''
     """
     samtools view -b -s ${sampling} ${chip_sample_id}.dedup.sorted.bam > ${chip_sample_id}.${sampling}.dedup.sorted.bam
     macs2 callpeak \\
         -t ${chip_sample_id}.${sampling}.dedup.sorted.bam \\
         $ctrl \\
         $broad \\
         -f BAM \\
         -g $REF_macs \\
         -n ${analysis_id}.${sampling} \\
         -q 0.01
     """
  }

  process saturation_r {
     tag "${saturation_results_collection[0].baseName}"
     publishDir "${params.outdir}/macs/saturation", mode: 'copy'

     input:
     file macsconfig from macsconfig
     file countstat from countstat_results
     file saturation_results_collection from saturation_results.collect()

     output:
     file '*.{txt,pdf}' into saturation_summary

     when: REF_macs

     script:
     """
     saturation_results_processing.r $params.rlocation $macsconfig $countstat $saturation_results_collection
     """
  }
}


/*
 * STEP 10 Post peak calling processing
 */

process chippeakanno {
    tag "${macs_peaks_collection[0].baseName}"
    publishDir "${params.outdir}/macs/chippeakanno", mode: 'copy'

    input:
    file macs_peaks_collection from macs_peaks.collect()
    file gtf from gtf

    output:
    file '*.{txt,bed}' into chippeakanno_results

    when: REF_macs

    script:
    filtering = params.blacklist_filtering ? "${params.blacklist}" : "No-filtering"
    """
    post_peak_calling_processing.r $params.rlocation $REF_macs $filtering $gtf $macs_peaks_collection
    """
}

/*
 * Parse software version numbers
 */
process get_software_versions {

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml

    script:
    """
    echo ${params.version} > v_ngi_chipseq.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    fastqc --version > v_fastqc.txt
    trim_galore --version > v_trim_galore.txt
    echo \$(bwa 2>&1) > v_bwa.txt
    samtools --version > v_samtools.txt
    bedtools --version > v_bedtools.txt
    echo "version" \$(java -Xmx2g -jar \$PICARD_HOME/picard.jar MarkDuplicates --version 2>&1) >v_picard.txt
    echo \$(plotFingerprint --version 2>&1) > v_deeptools.txt
    echo \$(ngs.plot.r 2>&1) > v_ngsplot.txt
    echo \$(macs2 --version 2>&1) > v_macs2.txt
    multiqc --version > v_multiqc.txt
    scrape_software_versions.py > software_versions_mqc.yaml
    """
}


/*
 * STEP 11 MultiQC
 */

process multiqc {
    tag "$prefix"
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file multiqc_config
    file (fastqc:'fastqc/*') from fastqc_results.collect()
    file ('trimgalore/*') from trimgalore_results.collect()
    file ('samtools/*') from samtools_stats.collect()
    file ('picard/*') from picard_reports.collect()
    file ('deeptools/*') from deepTools_multiqc.collect()
    file ('phantompeakqualtools/*') from spp_out_mqc.collect()
    file ('phantompeakqualtools/*') from calculateNSCRSC_results.collect()
    file ('software_versions/*') from software_versions_yaml.collect()

    output:
    file '*multiqc_report.html' into multiqc_report
    file '*_data' into multiqc_data
    file '.command.err' into multiqc_stderr
    val prefix into multiqc_prefix

    script:
    prefix = fastqc[0].toString() - '_fastqc.html' - 'fastqc/'
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    """
    multiqc -f $rtitle $rfilename --config $multiqc_config . 2>&1
    """
}

/*
 * STEP 12 - Output Description HTML
 */
process output_documentation {
    tag "$prefix"
    publishDir "${params.outdir}/Documentation", mode: 'copy'

    input:
    val prefix from multiqc_prefix
    file output from output_docs

    output:
    file "results_description.html"

    script:
    def rlocation = params.rlocation ?: ''
    """
    markdown_to_html.r $output results_description.html $rlocation
    """
}


/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/chipseq] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[nf-core/chipseq] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = params.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    if(workflow.container) email_fields['summary']['Docker image'] = workflow.container

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir" ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (params.email) {
        try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "[nf-core/chipseq] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[nf-core/chipseq] Sent summary e-mail to $params.email (mail)"
        }
    }

    // Switch the embedded MIME images with base64 encoded src
    nfcorechipseqlogo = new File("$baseDir/assets/nf-core_chipseq_logo.png").bytes.encodeBase64().toString()
    scilifelablogo = new File("$baseDir/assets/SciLifeLab_logo.png").bytes.encodeBase64().toString()
    ngilogo = new File("$baseDir/assets/NGI_logo.png").bytes.encodeBase64().toString()
    email_html = email_html.replaceAll(~/cid:nfcorechipseqlogo/, "data:image/png;base64,$nfcorechipseqlogo")
    email_html = email_html.replaceAll(~/cid:scilifelablogo/, "data:image/png;base64,$scilifelablogo")
    email_html = email_html.replaceAll(~/cid:ngilogo/, "data:image/png;base64,$ngilogo")

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/Documentation/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    log.info "[nf-core/chipseq] Pipeline Complete"

    if(!workflow.success){
        if( workflow.profile == 'standard'){
            if ( "hostname".execute().text.contains('.uppmax.uu.se') ) {
                log.error "====================================================\n" +
                        "  WARNING! You are running with the default 'standard'\n" +
                        "  pipeline config profile, which runs on the head node\n" +
                        "  and assumes all software is on the PATH.\n" +
                        "  This is probably why everything broke.\n" +
                        "  Please use `-profile uppmax` to run on UPPMAX clusters.\n" +
                        "============================================================"
            }
        }
    }
}
