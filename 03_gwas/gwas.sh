#!/usr/bin/env bash
# =============================================================================
# gwas.sh — Step 3: isoform-level GWAS with EMMAX (via gwasor.py)
# -----------------------------------------------------------------------------
# For each trait, gwasor.py runs the full `pri_gwas` chain:
#   phenotype keep -> quantile normalization -> genotype subset (maf 0.05)
#   -> EMMAX kinship (BN) -> EMMAX association -> .ps -> .qqman -> QQ/Manhattan
#
# Inputs :
#   $GT_PREFIX   — merged genotype prefix from Step 2 (needs .tfam/.map too)
#   $TRAITS      — space-separated trait names; each trait must have a
#                  phenotype file "<trait>" (two columns, tab: famID, value)
# Outputs:
#   gwasout_<trait>/*.qqman + *.ps + .man.jpeg/.qq.jpeg (per trait)
#
# Usage   : source ../config.env && bash gwas.sh
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
source "$REPO_ROOT/config.env"

GWAS_OUT="${GWAS_OUT:-./gwasout}"          # output directory (created per trait inside)
BFILE="${BFILE:-$GT_PREFIX}"               # genotype prefix (Step 2 output)
MAPFILE="${MAPFILE:-$DATA_DIR/map/genetic_map.map}"  # isoform genetic map (for .qqman conversion)

mkdir -p "$GWAS_OUT"

for T in $TRAITS; do
    echo ""
    echo "=== GWAS: $T ==="
    "$PYTHON" "$SCRIPT_DIR/gwasor.py" pri_gwas \
        -b "$BFILE" \
        -p "$T" \
        -t no \
        -sig "$SIG" \
        -o "$GWAS_OUT/gwasout_$T" \
        -map "$MAPFILE"
done
echo ""
echo "=== Step 3 done: results in $GWAS_OUT/gwasout_*/ (*.qqman) ==="
