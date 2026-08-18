#!/usr/bin/env bash
# =============================================================================
# plink.sh — Step 1b: convert the PAV matrix to a PLINK dataset
# -----------------------------------------------------------------------------
# 1) pav_xls2ped.py  : isoform PAV matrix (0/1) -> PLINK PED (AA/GG encoding),
#                      isoforms restricted to those present in the genetic map.
# 2) plink --make-bed: PED+MAP -> binary PLINK (bed/bim/fam)
# 3) plink recode12  : binary -> transpose format (tped/tfam) used by gwasor.py
#
# Inputs  : $PAV_MATRIX                 — isoform PAV matrix (samples as columns)
#           genetic_map.map            — produced by Step 1a (bed2map)
# Outputs : $PAV_PREFIX.{bed,bim,fam}  — binary PLINK dataset
#           $PAV_PREFIX.{tped,tfam}    — transpose-format genotype
#
# Usage   : source ../../config.env && bash plink.sh
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
source "$REPO_ROOT/config.env"

MAP_FILE="${MAP_FILE:-$DATA_DIR/map/genetic_map.map}"

# Step 1: PAV matrix -> PED (uses the genetic map to order/filter isoforms)
echo "=== Step 1b: PAV matrix -> PED ==="
"$PYTHON" "$SCRIPT_DIR/pav_xls2ped.py" "$PAV_MATRIX" "$MAP_FILE"

# PLINK needs a .map next to the .ped (same isoform order as the PED columns)
PED_BASE="$(basename "$PAV_MATRIX")"
cp "$MAP_FILE" "$PED_BASE.map"

echo "=== Step 1b: PED -> binary PLINK ==="
"$PLINK" --file "$PED_BASE" --chr-set 29 --make-bed --out "$PAV_PREFIX"

echo "=== Step 1b: binary -> transpose (tped/tfam) ==="
"$PLINK" --bfile "$PAV_PREFIX" --recode12 --output-missing-genotype 0 \
        --transpose --out "$PAV_PREFIX.maf0.05" --maf 0.05 --chr-set 29

echo "  $PAV_PREFIX.bed/.bim/.fam : $(wc -l < "$PAV_PREFIX".bim) variants, $(wc -l < "$PAV_PREFIX".fam) samples"
echo "  $PAV_PREFIX.maf0.05.tped/.tfam : transpose-format genotype"
echo "=== Step 1b done ==="
