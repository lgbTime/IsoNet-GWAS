#!/usr/bin/env bash
# =============================================================================
# get_uniq_pos.sh — Step 1a: build a unique-position isoform genetic map
# -----------------------------------------------------------------------------
# For every per-chromosome isoform BED file, deduplicate the isoform TSS
# positions so that each isoform gets a unique coordinate (required for a
# valid genetic map). Then concatenate all chromosomes into one map file.
#
# Inputs  : $BED_DIR/*.bed   (tab-separated, no header;
#                             col1=chr, col2=TSS, col3=TES, col4=isoform ID)
# Outputs : $MAP_OUT_DIR/out_*.bed        — per-chromosome unique-position BED
#           $MAP_OUT_DIR/genetic_map.bed  — concatenated, sorted genetic map
#           $MAP_OUT_DIR/genetic_map.map  — PLINK-ready map (chr, iso, 0, pos)
#
# Usage   : source ../config.env && bash get_uniq_pos.sh
# =============================================================================
set -euo pipefail

# Resolve repo root and load configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
source "$REPO_ROOT/config.env"

DEDUP_SCRIPT="${DEDUP_SCRIPT:-$SCRIPT_DIR/reassign_repeated_startpos2_v2.py}"
MAP_OUT_DIR="${MAP_OUT_DIR:-$DATA_DIR/map}"

mkdir -p "$MAP_OUT_DIR"
cd "$MAP_OUT_DIR"

if ! ls "$BED_DIR"/*.bed >/dev/null 2>&1; then
    echo "ERROR: no *.bed files found in BED_DIR=$BED_DIR" >&2
    exit 1
fi

echo "=== Step 1a: deduplicating isoform positions per chromosome ==="
for bed in "$BED_DIR"/*.bed; do
    echo "  -> $(basename "$bed")"
    "$PYTHON" "$DEDUP_SCRIPT" "$bed"
done

echo "=== Step 1a: concatenating chromosomes into genetic map ==="
cat out_*.bed | sort -k1,1 -k2,2n > genetic_map.bed

# PLINK .map format: chr  isoformID  genetic_cM  physical_pos
awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $4, 0, $2}' genetic_map.bed > genetic_map.map

echo "  genetic_map.bed: $(wc -l < genetic_map.bed) isoforms"
echo "  genetic_map.map : $(wc -l < genetic_map.map) markers"
echo "=== Step 1a done ==="
