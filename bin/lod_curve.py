#!/usr/bin/env python3
"""
lod_curve.py — SeraCare LOD sensitivity analysis and curve plot.

For each dilution point, checks how many of the 49 covered SeraCare positions
are detected by LoFreq and VarDict, computes sensitivity, fits a logistic curve,
and reports the LOD at 95% sensitivity.

Usage:
    python3 lod_curve.py \
        --vcf-dir    /path/to/vcf/output \
        --truth-set  seracare_v4_truth_set_hg38.csv \
        --samples    WT=SGULP16_9,2.5=SGULP16_1,1.25=SGULP16_2,0.625=SGULP16_3,0.3125=SGULP16_4 \
        --out-dir    /path/to/lod_results

--samples format: VAF_PCT=SAMPLE_NAME pairs, comma-separated.
  WT = negative control (used to flag artifacts).
  Other keys are expected VAF percentages (2.5, 1.25, 0.625, 0.3125).
"""

import argparse, os, sys, csv
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import norm as scipy_norm

CALLERS = ["lofreq", "vardict"]
PALETTE = {"lofreq": "#4472C4", "vardict": "#E87722"}
LOD_TARGET = 0.95   # 95% sensitivity threshold


# ── VCF parsing ───────────────────────────────────────────────────────────────

def parse_vcf(vcf_path):
    """Return list of dicts with chrom, pos (1-based), ref, alt, af."""
    records = []
    if not os.path.exists(vcf_path):
        return records
    with open(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, _, ref, alt = parts[0], int(parts[1]), parts[2], parts[3], parts[4]
            info = parts[7]
            af = None
            for field in info.split(";"):
                if field.startswith("AF="):
                    try:
                        af = float(field.split("=")[1].split(",")[0])
                    except ValueError:
                        pass
                    break
            if af is None and len(parts) >= 10:
                # VarDict puts AF in FORMAT/SAMPLE columns
                fmt = parts[8].split(":")
                smp = parts[9].split(":")
                if "AF" in fmt:
                    try:
                        af = float(smp[fmt.index("AF")])
                    except (ValueError, IndexError):
                        pass
            records.append({"chrom": chrom, "pos": pos, "ref": ref, "alt": alt, "af": af})
    return records


def called_set(vcf_records):
    """Return set of (chrom, pos) for quick lookup."""
    return {(r["chrom"], r["pos"]) for r in vcf_records}


# ── Logistic fit ──────────────────────────────────────────────────────────────

def logistic(x, k, x0):
    """Logistic function: sensitivity = 1 / (1 + exp(-k*(x - x0)))"""
    return 1.0 / (1.0 + np.exp(-k * (np.array(x) - x0)))


def fit_lod_curve(vafs, sensitivities):
    """
    Fit logistic curve; return (k, x0, lod_95).
    lod_95 = VAF at which sensitivity = 95%.
    Returns None values if fit fails or data is insufficient.
    """
    try:
        popt, _ = curve_fit(
            logistic, vafs, sensitivities,
            p0=[10, np.median(vafs)],
            bounds=([0, 0], [1000, max(vafs) * 10]),
            maxfev=10000
        )
        k, x0 = popt
        # Solve logistic for sensitivity = LOD_TARGET:
        # LOD_TARGET = 1 / (1 + exp(-k*(x - x0)))
        # => x = x0 - ln(1/LOD_TARGET - 1) / k
        if k > 0:
            lod_95 = x0 - np.log(1.0 / LOD_TARGET - 1.0) / k
        else:
            lod_95 = None
        return k, x0, lod_95
    except Exception:
        return None, None, None


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf-dir",   required=True)
    parser.add_argument("--truth-set", required=True,
                        help="CSV from seracare_v4_truth_set_hg38.csv")
    parser.add_argument("--samples",   required=True,
                        help="VAF_PCT=SAMPLE_NAME pairs, e.g. WT=S1,2.5=S2,...")
    parser.add_argument("--out-dir",   required=True)
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    # Parse sample map: {vaf_label -> sample_name}
    sample_map = {}
    for pair in args.samples.split(","):
        k, v = pair.strip().split("=")
        sample_map[k.strip()] = v.strip()

    expected_keys = {"WT", "2.5", "1.25", "0.625", "0.3125"}
    missing = expected_keys - set(sample_map.keys())
    if missing:
        print(f"WARNING: Missing sample keys: {missing}")

    dilutions_pct = [2.5, 1.25, 0.625, 0.3125]  # ordered high → low

    # Load truth set
    truth = pd.read_csv(args.truth_set)
    truth_positions = set(zip(truth["hg38_chrom"], truth["hg38_pos_vcf"]))
    N = len(truth_positions)
    print(f"Truth set: {N} positions")

    # ── Per-caller analysis ───────────────────────────────────────────────────
    results = {}   # {caller: {vaf_pct: {"detected": int, "sensitivity": float}}}
    wt_artifacts = {}  # {caller: set of (chrom,pos) called in WT}

    for caller in CALLERS:
        results[caller] = {}

        # WT negative control — identify artifacts
        wt_name = sample_map.get("WT")
        wt_calls = set()
        if wt_name:
            wt_vcf = os.path.join(args.vcf_dir, caller, f"{wt_name}_{caller}.vcf")
            wt_vcf_records = parse_vcf(wt_vcf)
            wt_calls = called_set(wt_vcf_records)
            wt_at_truth = wt_calls & truth_positions
            print(f"[{caller}] WT artifacts at truth positions: {len(wt_at_truth)}")
        wt_artifacts[caller] = wt_calls

        # Each dilution
        for vaf_pct in dilutions_pct:
            vaf_label = str(vaf_pct)
            sample_name = sample_map.get(vaf_label)
            if not sample_name:
                print(f"[{caller}] No sample for {vaf_pct}% — skipping")
                continue

            vcf_path = os.path.join(args.vcf_dir, caller, f"{sample_name}_{caller}.vcf")
            records = parse_vcf(vcf_path)
            calls = called_set(records)

            # True positives: in truth set, called in dilution, NOT in WT
            tp = (calls & truth_positions) - wt_calls
            sensitivity = len(tp) / N

            results[caller][vaf_pct] = {
                "sample":       sample_name,
                "vcf":          vcf_path,
                "total_calls":  len(calls),
                "tp":           len(tp),
                "fp_wt":        len(calls & wt_calls),  # calls that also appear in WT
                "sensitivity":  sensitivity,
            }
            print(f"[{caller}] {vaf_pct}%: {len(tp)}/{N} detected  (sensitivity={sensitivity:.1%})")

    # ── LOD curve fit ─────────────────────────────────────────────────────────
    print("\n── LOD fits ──")
    lod_summary = {}
    for caller in CALLERS:
        vafs = []
        sens = []
        for vaf_pct in dilutions_pct:
            if vaf_pct in results[caller]:
                vafs.append(vaf_pct)
                sens.append(results[caller][vaf_pct]["sensitivity"])
        if len(vafs) < 2:
            print(f"[{caller}] Not enough points for LOD fit")
            continue

        k, x0, lod_95 = fit_lod_curve(vafs, sens)
        lod_summary[caller] = {"k": k, "x0": x0, "lod_95": lod_95}
        if lod_95 is not None:
            print(f"[{caller}] LOD (95% sensitivity) = {lod_95:.3f}%")
        else:
            print(f"[{caller}] LOD fit failed")

    # ── Plot ──────────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle("SeraCare ctDNA Mutation Mix v4 — LOD Curve (MSK-ACCESS FLEX)", fontsize=13)

    # Left panel: sensitivity vs VAF (LOD curve)
    ax = axes[0]
    for caller in CALLERS:
        vafs = []
        sens = []
        for vaf_pct in dilutions_pct:
            if vaf_pct in results[caller]:
                vafs.append(vaf_pct)
                sens.append(results[caller][vaf_pct]["sensitivity"])
        if not vafs:
            continue

        color = PALETTE[caller]
        ax.scatter(vafs, [s * 100 for s in sens], color=color, s=80, zorder=5, label=f"{caller} (observed)")

        k, x0, lod_95 = lod_summary.get(caller, {}).get("k"), lod_summary.get(caller, {}).get("x0"), lod_summary.get(caller, {}).get("lod_95")
        if k is not None and x0 is not None:
            x_fit = np.linspace(0.1, max(vafs) * 1.5, 500)
            y_fit = logistic(x_fit, k, x0) * 100
            ax.plot(x_fit, y_fit, color=color, lw=2, label=f"{caller} (fit)")
            if lod_95 is not None and 0 < lod_95 < max(vafs) * 2:
                ax.axvline(lod_95, color=color, lw=1.2, ls="--", alpha=0.7)
                ax.text(lod_95 + 0.02, 10, f"{caller}\nLOD={lod_95:.2f}%",
                        color=color, fontsize=8, va="bottom")

    ax.axhline(95, color="grey", lw=1, ls=":", alpha=0.7, label="95% threshold")
    ax.set_xlabel("Expected VAF (%)")
    ax.set_ylabel("Sensitivity (%)")
    ax.set_ylim(0, 105)
    ax.set_title("LOD Curve — Sensitivity vs VAF")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # Right panel: per-variant detection heatmap at each dilution
    ax2 = axes[1]
    truth_sorted = truth.sort_values(["gene", "mutation"]).reset_index(drop=True)
    labels = [f"{r.gene} {r.mutation}" for _, r in truth_sorted.iterrows()]

    caller = CALLERS[0]  # show LoFreq on right panel
    heatmap_data = []
    col_labels = []
    for vaf_pct in dilutions_pct:
        sample_name = sample_map.get(str(vaf_pct))
        if not sample_name:
            continue
        vcf_path = os.path.join(args.vcf_dir, caller, f"{sample_name}_{caller}.vcf")
        calls = called_set(parse_vcf(vcf_path)) - wt_artifacts[caller]
        row = []
        for _, r in truth_sorted.iterrows():
            detected = 1 if (r["hg38_chrom"], r["hg38_pos_vcf"]) in calls else 0
            row.append(detected)
        heatmap_data.append(row)
        col_labels.append(f"{vaf_pct}%")

    if heatmap_data:
        mat = np.array(heatmap_data)  # shape: (n_dilutions, n_variants)
        im = ax2.imshow(mat, aspect="auto", cmap="RdYlGn", vmin=0, vmax=1, interpolation="none")
        ax2.set_yticks(range(len(col_labels)))
        ax2.set_yticklabels(col_labels, fontsize=9)
        ax2.set_xticks(range(len(labels)))
        ax2.set_xticklabels(labels, rotation=90, fontsize=6)
        ax2.set_title(f"Per-variant detection — {caller} (green=detected, red=missed)")
        plt.colorbar(im, ax=ax2, shrink=0.6)

    plt.tight_layout()
    out_png = os.path.join(args.out_dir, "lod_curve.png")
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    print(f"\nPlot saved: {out_png}")

    # ── Summary CSV ──────────────────────────────────────────────────────────
    rows = []
    for caller in CALLERS:
        for vaf_pct in dilutions_pct:
            if vaf_pct not in results[caller]:
                continue
            r = results[caller][vaf_pct]
            lod_95 = lod_summary.get(caller, {}).get("lod_95")
            rows.append({
                "caller": caller,
                "expected_vaf_pct": vaf_pct,
                "sample": r["sample"],
                "variants_detected": r["tp"],
                "variants_total": N,
                "sensitivity_pct": round(r["sensitivity"] * 100, 1),
                "total_vcf_calls": r["total_calls"],
                "wt_artifact_overlap": r["fp_wt"],
                "lod_95_pct": round(lod_95, 3) if lod_95 else "NA",
            })

    summary_csv = os.path.join(args.out_dir, "lod_summary.csv")
    if rows:
        pd.DataFrame(rows).to_csv(summary_csv, index=False)
        print(f"Summary CSV: {summary_csv}")
        print("\n" + pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    main()
