#!/usr/bin/env bash
# =============================================================================
# gwas_norm.sh — Step 3 variant: GWAS with a different output prefix
# -----------------------------------------------------------------------------
# Identical to gwas.sh but writes to norm_gwasout_<trait>/, e.g. for running
# with a normalized phenotype transformation.
#
# Usage   : source ../config.env && bash gwas_norm.sh
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
source "$REPO_ROOT/config.env"

GWAS_OUT="${GWAS_OUT:-./gwasout_norm}"     # output directory (created per trait inside)
BFILE="${BFILE:-$GT_PREFIX}"
MAPFILE="${MAPFILE:-$DATA_DIR/map/genetic_map.map}"

mkdir -p "$GWAS_OUT"

for T in $TRAITS; do
    echo ""
    echo "=== GWAS (norm): $T ==="
    "$PYTHON" "$SCRIPT_DIR/gwasor.py" pri_gwas \
        -b "$BFILE" \
        -p "$T" \
        -t no \
        -sig "$SIG" \
        -o "$GWAS_OUT/norm_gwasout_$T" \
        -map "$MAPFILE"
done
echo ""
echo "=== Step 3 (norm) done: results in $GWAS_OUT/norm_gwasout_*/ ==="
