#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ── Parameters ────────────────────────────────────────────────────────────────
params.bam_dir  = null   // required: S3 folder containing *_aligned_sorted.bam files
params.ref      = null   // required: S3 path to hg38 FASTA (+ .fai in same folder)
params.bed      = 's3://volta-labs-sequencing-data/reference-genome-index/hybrid-capture-target-files/SG_MSKACCESSFLEX_v1_target_hg38.bed'
params.outdir   = null   // required: S3 output path
params.min_vaf  = 0.001  // minimum VAF threshold (0.1%)
params.min_cov  = 5      // minimum read depth to call

if (!params.bam_dir) { error "Please provide --bam_dir <s3://bucket/path/to/bams/>" }
if (!params.ref)     { error "Please provide --ref <s3://bucket/path/to/hg38.fasta>" }
if (!params.outdir)  { error "Please provide --outdir <s3://bucket/path/>" }

// ── Process: LoFreq ───────────────────────────────────────────────────────────
process LOFREQ {
    tag "$sample"

    conda 'bioconda::lofreq bioconda::samtools'

    publishDir "${params.outdir}/lofreq/", mode: 'copy', overwrite: true

    cpus   8
    memory { 16.GB * task.attempt }
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3

    input:
    tuple val(sample), path(bam), path(bai)
    path(ref)
    path(bed)

    output:
    tuple val(sample), path("${sample}_lofreq.vcf"), emit: vcf

    script:
    def ref_file = ref.find { it.name.endsWith('.fasta') || it.name.endsWith('.fa') }
    """
    # Normalise BED chromosome names (add chr prefix if absent)
    awk 'BEGIN{OFS="\\t"} /^[^#]/ && \$2 !~ /^[0-9]/ {next} /^[^#]/ && \$1 !~ /^chr/ {\$1="chr"\$1} {print}' \\
        ${bed} > normalised.bed

    # Recalibrate indel qualities
    lofreq indelqual --dindel -f ${ref_file} ${bam} -o indelqual.bam
    samtools index indelqual.bam

    # Call variants (single-threaded; call-parallel has a macOS bgzip bug)
    lofreq call \\
        -f ${ref_file} \\
        -l normalised.bed \\
        --call-indels \\
        --min-cov ${params.min_cov} \\
        -o ${sample}_lofreq.vcf \\
        indelqual.bam

    rm -f indelqual.bam indelqual.bam.bai
    """
}

// ── Process: VarDict ─────────────────────────────────────────────────────────
process VARDICT {
    tag "$sample"

    conda 'bioconda::vardict-java'

    publishDir "${params.outdir}/vardict/", mode: 'copy', overwrite: true

    cpus   8
    memory { 16.GB * task.attempt }
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3

    input:
    tuple val(sample), path(bam), path(bai)
    path(ref)
    path(bed)

    output:
    tuple val(sample), path("${sample}_vardict.vcf"), emit: vcf

    script:
    def ref_file = ref.find { it.name.endsWith('.fasta') || it.name.endsWith('.fa') }
    """
    # Normalise BED chromosome names
    awk 'BEGIN{OFS="\\t"} /^[^#]/ && \$2 !~ /^[0-9]/ {next} /^[^#]/ && \$1 !~ /^chr/ {\$1="chr"\$1} {print}' \\
        ${bed} > normalised.bed

    vardict-java \\
        -G ${ref_file} \\
        -f ${params.min_vaf} \\
        -N ${sample} \\
        -b ${bam} \\
        -c 1 -S 2 -E 3 -g 4 \\
        -th ${task.cpus} \\
        normalised.bed \\
    | teststrandbias.R \\
    | var2vcf_valid.pl \\
        -N ${sample} \\
        -E \\
        -f ${params.min_vaf} \\
    > ${sample}_vardict.vcf
    """
}

// ── Workflow ──────────────────────────────────────────────────────────────────
workflow {

    // Auto-discover BAMs from bam_dir
    Channel
        .fromPath("${params.bam_dir}/*_aligned_sorted.bam")
        .map { bam ->
            def sample = bam.name.replaceAll(/_aligned_sorted\.bam$/, '')
            def bai    = bam.parent.resolve("${bam.name}.bai")
            tuple(sample, bam, bai)
        }
        .set { bams_ch }

    // Reference files (FASTA + .fai)
    ref_ch = Channel.fromPath("${params.ref}*").collect()

    // BED file
    bed_ch = Channel.fromPath(params.bed)

    // Run both callers in parallel for every sample
    LOFREQ(bams_ch, ref_ch, bed_ch)
    VARDICT(bams_ch, ref_ch, bed_ch)
}
