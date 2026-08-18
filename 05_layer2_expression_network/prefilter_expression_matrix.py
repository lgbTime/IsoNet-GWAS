#!/usr/bin/env python3
"""
prefilter_expression_matrix.py — Filter isoform-level TPM matrix to remove
                                  low-expression noise before downstream analysis.

Rationale:
    The raw isoform TPM matrix (~712k isoforms × 240 samples) is dominated by
    noise: 50% of isoforms have mean TPM < 0.01 and only ~7% are confidently
    expressed (TPM ≥ 1 in ≥ 20% of samples).  Running co-expression or TF-target
    correlation on unfiltered data inflates the noise floor, produces spurious
    edges, and dramatically increases I/O and compute cost.

    This script applies isoform-level filters *before* gene-level aggregation,
    matching the standard used in Nature/Science plant co-expression studies.

Filters (isoform passes if ANY condition is met — UNION logic):
    1. Mean TPM across all samples ≥ --min-mean-tpm     (default 0.5)
    2. TPM ≥ --expr-threshold in ≥ --min-samples-pct of samples   (default 1.0, 20%)
    3. (Optional) TPM ≥ --expr-threshold in ≥ --min-samples-abs samples

Output:
    A cleaned TPM file with the same format (geneID:transcriptID<TAB>sample1<TAB>...)
    containing only isoforms that pass at least one filter.

Usage:
    # Default filter (mean TPM ≥ 0.5 OR TPM ≥ 1 in ≥ 20% of samples)
    python3 prefilter_expression_matrix.py \
        --input  Pop_isoform_TPM_by_compareGTF.txt \
        --output Pop_isoform_TPM_filtered.txt

    # Stringent filter (mean TPM ≥ 1.0 AND TPM ≥ 1 in ≥ 50% of samples)
    python3 prefilter_expression_matrix.py \
        --input  Pop_isoform_TPM_by_compareGTF.txt \
        --output Pop_isoform_TPM_stringent.txt \
        --min-mean-tpm 1.0 --min-samples-pct 50 --mode intersection

    # Keep tissue-specific isoforms (low mean TPM but high in ≥ 10 samples)
    python3 prefilter_expression_matrix.py \
        --input  Pop_isoform_TPM_by_compareGTF.txt \
        --output Pop_isoform_TPM_tissuespecific.txt \
        --min-mean-tpm 0.1 --expr-threshold 2.0 --min-samples-abs 10
"""

import argparse
import os
import sys
import time
import numpy as np


def parse_args():
    p = argparse.ArgumentParser(
        description="Filter isoform-level TPM matrix by expression level"
    )
    p.add_argument("--input", "-i", required=True,
                   help="Input TPM matrix (geneID:transcriptID <TAB> sample1 sample2 ...)")
    p.add_argument("--output", "-o", required=True,
                   help="Output filtered TPM matrix (same format)")
    p.add_argument("--min-mean-tpm", dest="min_mean_tpm",
                   type=float, default=0.5,
                   help="Minimum mean TPM across all samples [default: 0.5]")
    p.add_argument("--expr-threshold", dest="expr_threshold",
                   type=float, default=1.0,
                   help="TPM value that counts as 'expressed' for sample-count filter [default: 1.0]")
    p.add_argument("--min-samples-pct", dest="min_samples_pct",
                   type=float, default=20.0,
                   help="Minimum %% of samples with TPM ≥ expr-threshold [default: 20]")
    p.add_argument("--min-samples-abs", dest="min_samples_abs",
                   type=int, default=0,
                   help="Minimum absolute number of samples with TPM ≥ expr-threshold "
                        "(overrides --min-samples-pct if > 0)")
    p.add_argument("--mode", choices=["union", "intersection"], default="union",
                   help="How to combine filters.  'union' = pass if ANY filter passes "
                        "(permissive, keeps tissue-specific genes).  'intersection' = "
                        "pass only if ALL filters pass (stringent) [default: union]")
    p.add_argument("--quiet", "-q", action="store_true",
                   help="Suppress progress output")
    return p.parse_args()


def main():
    args = parse_args()

    t_start = time.time()

    # ---- Read header ----
    with open(args.input, "r") as fh:
        header = fh.readline().rstrip("\n")
    sample_names = header.split("\t")[1:]
    n_samples = len(sample_names)

    if not args.quiet:
        print(f"Input:  {args.input}")
        print(f"        {n_samples} samples")
        print(f"Filters:")
        print(f"  (A) mean TPM ≥ {args.min_mean_tpm}")
        if args.min_samples_abs > 0:
            print(f"  (B) TPM ≥ {args.expr_threshold} in ≥ {args.min_samples_abs} samples "
                  f"({100*args.min_samples_abs/n_samples:.1f}%)")
        else:
            min_n = max(1, int(round(n_samples * args.min_samples_pct / 100.0)))
            print(f"  (B) TPM ≥ {args.expr_threshold} in ≥ {min_n} samples "
                  f"({args.min_samples_pct:.0f}%)")
        print(f"  Mode: {args.mode.upper()}")
        print(f"Reading & filtering …")

    # ---- Determine sample-count threshold ----
    if args.min_samples_abs > 0:
        min_samples_n = args.min_samples_abs
    else:
        min_samples_n = max(1, int(round(n_samples * args.min_samples_pct / 100.0)))

    # ---- Single-pass filter ----
    n_total = 0
    n_pass_mean = 0
    n_pass_count = 0
    n_pass_both = 0
    n_pass_any = 0

    kept_lines = [header]

    with open(args.input, "r") as fh:
        fh.readline()  # skip header, already read
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue

            n_total += 1
            isoform_id = parts[0]

            try:
                values = np.array([float(v) for v in parts[1:]], dtype=np.float32)
            except ValueError:
                # Unparseable row — keep it (could be metadata or annotation row)
                kept_lines.append(line)
                continue

            # Filter A: mean TPM
            mean_tpm = float(values.mean())
            pass_mean = mean_tpm >= args.min_mean_tpm
            if pass_mean:
                n_pass_mean += 1

            # Filter B: sample count above expression threshold
            n_expr_samples = int(np.sum(values >= args.expr_threshold))
            pass_count = n_expr_samples >= min_samples_n
            if pass_count:
                n_pass_count += 1

            # Combine
            if pass_mean and pass_count:
                n_pass_both += 1

            if args.mode == "union":
                keep = pass_mean or pass_count
            else:  # intersection
                keep = pass_mean and pass_count

            if keep:
                n_pass_any += 1
                kept_lines.append(line)

    # ---- Write output ----
    with open(args.output, "w") as fh:
        for l in kept_lines:
            fh.write(l + "\n")

    # ---- Report ----
    elapsed = time.time() - t_start
    if not args.quiet:
        pct_pass = 100.0 * n_pass_any / n_total if n_total > 0 else 0.0
        print(f"\n{'='*60}")
        print(f"Filtering complete ({elapsed:.1f}s)")
        print(f"{'='*60}")
        print(f"  Total isoforms in:        {n_total:>10,}")
        print(f"  Passed mean TPM ≥ {args.min_mean_tpm}:     {n_pass_mean:>10,}  "
              f"({100*n_pass_mean/max(n_total,1):.1f}%)")
        print(f"  Passed sample-count:      {n_pass_count:>10,}  "
              f"({100*n_pass_count/max(n_total,1):.1f}%)")
        print(f"  Passed BOTH filters:      {n_pass_both:>10,}  "
              f"({100*n_pass_both/max(n_total,1):.1f}%)")
        print(f"  ─────────────────────────────────────────")
        print(f"  Retained ({args.mode}):        {n_pass_any:>10,}  "
              f"({pct_pass:.1f}%)")
        print(f"  Removed:                  {n_total - n_pass_any:>10,}  "
              f"({100*(n_total-n_pass_any)/max(n_total,1):.1f}%)")
        print(f"\nOutput: {args.output}")
        print(f"        {n_pass_any:,} isoforms × {n_samples} samples")
        print(f"        {os.path.getsize(args.output)/1024/1024:.1f} MB"
              if "os" in dir() else "")

    return 0


if __name__ == "__main__":
    sys.exit(main())
