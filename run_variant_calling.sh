#!/usr/bin/env bash
# run_variant_calling.sh — Run LoFreq and VarDict on SeraCare LOD BAMs
#
# Usage:
#   bash run_variant_calling.sh \
#       --bam-dir  /path/to/bam/folder \
#       --ref      /path/to/hg38.fasta \
#       --bed      /path/to/SG_MSKACCESSFLEX_v1_target_hg38.bed \
#       --out-dir  /path/to/output
#
# Input BAMs must follow illumina-metrics-nf naming: {sample}_aligned_sorted.bam
# Requires: lofreq, vardict-java, samtools, bcftools (all in lod-pipeline conda env)
#
# Expected samples (edit SAMPLES array to match your run):
#   WT        = negative control (SeraCare WT)
#   2p5pct    = 2.5% dilution
#   1p25pct   = 1.25% dilution
#   0p625pct  = 0.625% dilution
#   0p3125pct = 0.3125% dilution
# -----------------------------------------------------------------------------
set -euo pipefail

# ── Parse arguments ───────────────────────────────────────────────────────────
BAM_DIR=""
REF=""
BED=""
OUT_DIR=""
MIN_VAF=0.001   # 0.1% — catch even the lowest dilution point
MIN_AF_VAR=0.001
THREADS=4

while [[ $# -gt 0 ]]; do
    case $1 in
        --bam-dir)  BAM_DIR="$2"; shift 2 ;;
        --ref)      REF="$2";     shift 2 ;;
        --bed)      BED="$2";     shift 2 ;;
        --out-dir)  OUT_DIR="$2"; shift 2 ;;
        --threads)  THREADS="$2"; shift 2 ;;
        *) echo "Unknown arg: $1"; exit 1 ;;
    esac
done

if [[ -z "$BAM_DIR" || -z "$REF" || -z "$BED" || -z "$OUT_DIR" ]]; then
    echo "ERROR: --bam-dir, --ref, --bed, and --out-dir are all required."
    exit 1
fi

mkdir -p "$OUT_DIR/lofreq" "$OUT_DIR/vardict" "$OUT_DIR/logs"

echo "=== Variant Calling — SeraCare LOD ==="
echo "BAM dir : $BAM_DIR"
echo "Ref     : $REF"
echo "BED     : $BED"
echo "Out dir : $OUT_DIR"
echo "Min VAF : $MIN_VAF"
echo ""

# ── Discover BAMs ─────────────────────────────────────────────────────────────
BAMS=("$BAM_DIR"/*_aligned_sorted.bam)
if [[ ${#BAMS[@]} -eq 0 || ! -f "${BAMS[0]}" ]]; then
    echo "ERROR: No *_aligned_sorted.bam files found in $BAM_DIR"
    exit 1
fi

echo "Found ${#BAMS[@]} BAM(s):"
for b in "${BAMS[@]}"; do echo "  $(basename $b)"; done
echo ""

# ── Normalise BED: ensure chr prefix ─────────────────────────────────────────
NORM_BED="$OUT_DIR/normalised_design.bed"
awk 'BEGIN{OFS="\t"} /^[^#]/ && $2 !~ /^[0-9]/ {next} /^[^#]/ && $1 !~ /^chr/ {$1="chr"$1} {print}' \
    "$BED" > "$NORM_BED"
echo "Normalised BED: $(wc -l < "$NORM_BED") intervals"
echo ""

# ── Process each BAM ─────────────────────────────────────────────────────────
for BAM in "${BAMS[@]}"; do
    SAMPLE=$(basename "$BAM" _aligned_sorted.bam)
    echo "──────────────────────────────────────"
    echo "Sample: $SAMPLE"
    LOG="$OUT_DIR/logs/${SAMPLE}.log"

    # Index BAM if needed
    if [[ ! -f "${BAM}.bai" ]]; then
        echo "  Indexing BAM..."
        samtools index -@ "$THREADS" "$BAM" 2>>"$LOG"
    fi

    # ── LoFreq ────────────────────────────────────────────────────────────────
    LF_OUT="$OUT_DIR/lofreq/${SAMPLE}_lofreq.vcf"
    if [[ -f "$LF_OUT" ]]; then
        echo "  [LoFreq] Already done, skipping."
    else
        echo "  [LoFreq] Running indelqual recalibration..."
        INDELQUAL_BAM="$OUT_DIR/lofreq/${SAMPLE}_indelqual.bam"
        lofreq indelqual --dindel \
            -f "$REF" \
            "$BAM" \
            -o "$INDELQUAL_BAM" \
            >> "$LOG" 2>&1

        samtools index "$INDELQUAL_BAM" >> "$LOG" 2>&1

        echo "  [LoFreq] Calling variants (single-threaded, ~5-10 min)..."
        # Use lofreq call directly — call-parallel has a macOS bgzip bug in its
        # internal filter step that causes it to fail regardless of flags.
        lofreq call \
            -f "$REF" \
            -l "$NORM_BED" \
            --call-indels \
            --min-cov 5 \
            -o "$LF_OUT" \
            "$INDELQUAL_BAM" \
            >> "$LOG" 2>&1

        rm -f "$INDELQUAL_BAM" "${INDELQUAL_BAM}.bai"
        N=$(grep -cv "^#" "$LF_OUT" 2>/dev/null || echo 0)
        echo "  [LoFreq] Done — $N variants called"
    fi

    # ── VarDict ───────────────────────────────────────────────────────────────
    VD_RAW="$OUT_DIR/vardict/${SAMPLE}_vardict_raw.txt"
    VD_OUT="$OUT_DIR/vardict/${SAMPLE}_vardict.vcf"

    if [[ -f "$VD_OUT" ]]; then
        echo "  [VarDict] Already done, skipping."
    else
        echo "  [VarDict] Calling variants..."
        vardict-java \
            -G "$REF" \
            -f "$MIN_AF_VAR" \
            -N "$SAMPLE" \
            -b "$BAM" \
            -c 1 -S 2 -E 3 -g 4 \
            -th "$THREADS" \
            "$NORM_BED" \
            > "$VD_RAW" 2>>"$LOG"

        # VarDict post-processing pipeline
        cat "$VD_RAW" \
            | teststrandbias.R \
            | var2vcf_valid.pl \
                -N "$SAMPLE" \
                -E \
                -f "$MIN_AF_VAR" \
            > "$VD_OUT" 2>>"$LOG"

        N=$(grep -cv "^#" "$VD_OUT" 2>/dev/null || echo 0)
        echo "  [VarDict] Done — $N variants called"
    fi

    echo ""
done

echo "======================================"
echo "All samples complete."
echo "Outputs:"
echo "  LoFreq VCFs : $OUT_DIR/lofreq/"
echo "  VarDict VCFs: $OUT_DIR/vardict/"
echo ""
echo "Next step: run lod_curve.py"
echo "  python3 lod_curve.py \\"
echo "      --vcf-dir    $OUT_DIR \\"
echo "      --truth-set  seracare_v4_truth_set_hg38.csv \\"
echo "      --samples    WT=WT_SAMPLE_NAME,2.5=SAMPLE_2P5,1.25=SAMPLE_1P25,0.625=SAMPLE_0P625,0.3125=SAMPLE_0P3125 \\"
echo "      --out-dir    $OUT_DIR/lod_results"
