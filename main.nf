#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ── Parameters ────────────────────────────────────────────────────────────────
params.bam_dir     = null   // required: S3 folder containing *.bam and *.bam.bai files
params.ref         = null   // required: S3 path to hg38 FASTA (+ .fai in same folder)
params.bed         = null   // required: S3 path to target BED (hg38)
params.gc_cache    = null   // required: S3 path to GC content cache CSV (chrom,start,end,gc)
params.samplesheet = null   // required: S3 path to samplesheet CSV (sample,condition,pool)
params.truth_set   = null   // optional: S3 path to truth set CSV; defaults to bundled file
params.outdir      = null   // required: S3 output path
params.min_vaf     = 0.001  // minimum VAF threshold (0.1%)
params.min_cov     = 5      // minimum read depth to call

if (!params.bam_dir)     { error "Please provide --bam_dir <s3://bucket/path/to/bams/>" }
if (!params.ref)         { error "Please provide --ref <s3://bucket/path/to/hg38.fasta>" }
if (!params.bed)         { error "Please provide --bed <s3://bucket/path/to/target.bed>" }
if (!params.gc_cache)    { error "Please provide --gc_cache <s3://bucket/path/to/gc_cache.csv>" }
if (!params.samplesheet) { error "Please provide --samplesheet <s3://bucket/path/to/samplesheet.csv>" }
if (!params.outdir)      { error "Please provide --outdir <s3://bucket/path/>" }

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

    # Call variants in parallel (Linux/AWS Batch only — call-parallel has a bgzip bug on macOS)
    lofreq call-parallel \\
        --pp-threads ${task.cpus} \\
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

// ── Process: LOD Curve ────────────────────────────────────────────────────────
process LOD_CURVE {
    tag "lod_analysis"

    conda 'conda-forge::python>=3.10 conda-forge::pandas conda-forge::matplotlib conda-forge::scipy conda-forge::numpy'

    publishDir "${params.outdir}/lod_results/", mode: 'copy', overwrite: true

    cpus   2
    memory '8 GB'
    time   '30 min'

    input:
    path(lofreq_vcfs, stageAs: "staged_lofreq/*")
    path(vardict_vcfs, stageAs: "staged_vardict/*")
    path(samplesheet)
    path(truth_set)
    path(lod_script)   // staged explicitly from bin/ — avoids PATH issues in Wave containers

    output:
    path("pool_*/lod_curve.png")
    path("pool_*/lod_summary.csv")

    script:
    """
    # Build vcf-dir structure expected by lod_curve.py
    mkdir -p vcf_dir/lofreq vcf_dir/vardict
    cp staged_lofreq/*.vcf vcf_dir/lofreq/ 2>/dev/null || true
    cp staged_vardict/*.vcf vcf_dir/vardict/ 2>/dev/null || true

    # Run LOD curve analysis once per pool
    python3 << 'PYEOF'
import csv, os, subprocess, sys

pools = {}
with open("${samplesheet}") as f:
    for row in csv.DictReader(f):
        if row["condition"].strip().lower() == "ignore":
            continue
        pool = row.get("pool", "1").strip()
        if pool not in pools:
            pools[pool] = []
        pools[pool].append(f"{row['condition'].strip()}={row['sample'].strip()}")

for pool, pairs in sorted(pools.items()):
    out_dir = f"pool_{pool}"
    os.makedirs(out_dir, exist_ok=True)
    samples_str = ",".join(pairs)
    print(f"Pool {pool}: {samples_str}", flush=True)
    subprocess.run([
        sys.executable, "${lod_script}",
        "--vcf-dir",   "vcf_dir",
        "--truth-set", "${truth_set}",
        "--samples",   samples_str,
        "--out-dir",   out_dir,
    ], check=True)
PYEOF
    """
}

// ── Workflow ──────────────────────────────────────────────────────────────────
workflow {

    // Load samplesheet into a set of sample names to process
    def sample_set = [] as Set
    file(params.samplesheet).readLines().drop(1).each { line ->
        def cols = line.split(',')
        if (cols.size() >= 1) sample_set << cols[0].trim()
    }

    // Auto-discover BAMs; filter to only samples in samplesheet
    // Strip _aligned_sorted.bam suffix and _S\d+ Illumina index suffix (matches coverage-plots-nf)
    Channel
        .fromPath("${params.bam_dir}/*.bam")
        .map { bam ->
            def sample = bam.name
                .replaceAll(/_aligned_sorted\.bam$/, '')
                .replaceAll(/_S\d+$/, '')
            def bai = bam.parent.resolve("${bam.name}.bai")
            tuple(sample, bam, bai)
        }
        .filter { sample, bam, bai -> sample_set.contains(sample) }
        .set { bams_ch }

    // Reference files (FASTA + .fai)
    ref_ch = Channel.fromPath("${params.ref}*").collect()

    // BED and GC cache
    bed_ch = file(params.bed)

    // Samplesheet and truth set
    samplesheet_ch = file(params.samplesheet)
    truth_set_ch   = params.truth_set
        ? file(params.truth_set)
        : file("${projectDir}/seracare_v4_truth_set_hg38.csv")

    // LOD curve script — staged explicitly so PATH is not needed in Wave containers
    lod_script_ch = file("${projectDir}/bin/lod_curve.py")

    // Run both callers in parallel for every sample
    LOFREQ(bams_ch, ref_ch, bed_ch)
    VARDICT(bams_ch, ref_ch, bed_ch)

    // Collect all VCFs then run LOD curve analysis
    lofreq_vcfs_ch = LOFREQ.out.vcf.map { sample, vcf -> vcf }.collect()
    vardict_vcfs_ch = VARDICT.out.vcf.map { sample, vcf -> vcf }.collect()

    LOD_CURVE(lofreq_vcfs_ch, vardict_vcfs_ch, samplesheet_ch, truth_set_ch, lod_script_ch)
}
