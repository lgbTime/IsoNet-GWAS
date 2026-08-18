#!/usr/bin/env bash
# =============================================================================
# run.sh — Step 2: merge extra marker VCFs (gene-body SVs + GLK PAV) into the
#          isoform-PAV PLINK dataset produced by Step 1.
# -----------------------------------------------------------------------------
# Inputs:
#   $GENE_BODY_VCF  — gene-body SV variants (e.g. extracted.gene_body.filtered_num_chr.vcf)
#   $GLK_PAV_VCF    — GLK PAV variant VCF (e.g. GLK.PAV.txt.GT.vcf)
#   $PAV_PREFIX     — PLINK dataset from Step 1b (lettome.bed/.bim/.fam)
#
# Only samples overlapping between the VCFs and the PAV dataset are kept;
# non-overlapping samples are dropped.
#
# Output: $GT_PREFIX.{bed,bim,fam,tped,tfam,map}  (binary + text PLINK)
#
# GT conversion (PLINK standard, consistent with add_new_pos2tpedFormat.py):
#   VCF diploid:  0/0->"1 1"  0/1->"1 2"  1/1->"2 2"  ./.->"0 0"
#   VCF haploid:  0->"1 1"    1->"2 2"    .->"0 0"
#   Chromosome:   Lsat_Salinas_v11_chr4 -> 4
#   Variant IDs:  variantID:GENE_ID   (e.g. 1_42099:LOC111905264)
#
# Usage   : source ../config.env && bash run.sh
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
source "$REPO_ROOT/config.env"

if ! command -v "$PLINK" >/dev/null 2>&1; then
    echo "ERROR: PLINK binary not found: $PLINK" >&2
    exit 1
fi

# -------------------------------------------------------------------
# Step 1: Convert main gene-body VCF -> PLINK tped/map/tfam
# -------------------------------------------------------------------
echo "=== Step 2.1: Converting gene-body VCF to PLINK ==="
"$PYTHON" "$SCRIPT_DIR/vcf_sv_to_plink.py" \
    --vcf  "$GENE_BODY_VCF" \
    --tfam "$PAV_PREFIX.tfam" \
    --out  vcf_variants

echo "=== Step 2.1b: VCF variants -> binary ==="
"$PLINK" \
    --tped vcf_variants.tped \
    --tfam vcf_variants.tfam \
    --make-bed \
    --out vcf_variants \
    --allow-no-sex \
    --chr-set 50 \
    --allow-extra-chr

# -------------------------------------------------------------------
# Step 2: Convert GLK PAV VCF -> PLINK tped/map/tfam
# -------------------------------------------------------------------
echo "=== Step 2.2: Converting GLK.PAV VCF to PLINK ==="
"$PYTHON" "$SCRIPT_DIR/vcf_sv_to_plink.py" \
    --vcf  "$GLK_PAV_VCF" \
    --tfam "$PAV_PREFIX.tfam" \
    --out  glk_vcf_variants

echo "=== Step 2.2b: GLK variants -> binary ==="
"$PLINK" \
    --tped glk_vcf_variants.tped \
    --tfam glk_vcf_variants.tfam \
    --make-bed \
    --out glk_vcf_variants \
    --allow-no-sex \
    --chr-set 50 \
    --allow-extra-chr

# -------------------------------------------------------------------
# Step 3: Subset existing PAV dataset to overlap samples
# -------------------------------------------------------------------
echo "=== Step 2.3: Subsetting $PAV_PREFIX to overlap samples ==="
"$PLINK" \
    --bfile "$PAV_PREFIX" \
    --keep vcf_variants.overlap_samples.txt \
    --make-bed \
    --out pav_overlap \
    --allow-no-sex \
    --chr-set 50 \
    --allow-extra-chr

echo "pav_overlap: $(wc -l < pav_overlap.bim) variants, $(wc -l < pav_overlap.fam) samples"

# -------------------------------------------------------------------
# Step 4: Merge all three — PAV + gene-body VCF + GLK VCF
# -------------------------------------------------------------------
echo "=== Step 2.4: Merging PAV + gene-body VCF ==="
"$PLINK" \
    --bfile pav_overlap \
    --bmerge vcf_variants \
    --make-bed \
    --out tmp_merged \
    --allow-no-sex \
    --chr-set 50 \
    --allow-extra-chr

echo "=== Step 2.4b: Merging tmp_merged + GLK VCF ==="
"$PLINK" \
    --bfile tmp_merged \
    --bmerge glk_vcf_variants \
    --make-bed \
    --out "$GT_PREFIX" \
    --allow-no-sex \
    --chr-set 50 \
    --allow-extra-chr

# -------------------------------------------------------------------
# Step 5: Re-export as tped/tfam/map so text-format files are in sync
# -------------------------------------------------------------------
echo "=== Step 2.5: Re-exporting tped/tfam from binary ==="
"$PLINK" \
    --bfile "$GT_PREFIX" \
    --recode transpose \
    --out "$GT_PREFIX" \
    --allow-no-sex \
    --chr-set 50 \
    --allow-extra-chr

# Generate .map from .bim (drop allele columns, tab-delimited)
awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $2, $3, $4}' \
    "$GT_PREFIX.bim" > "$GT_PREFIX.map"

# Make matching .fam from .tfam (PLINK recode writes .tfam space-delimited, but
# PLINK --make-bed writes .fam space-delimited — both are correct PLINK standard)
cp "$GT_PREFIX.tfam" "$GT_PREFIX.fam"

echo ""
echo "=== Step 2 done ==="
echo "Final files (space-delimited .fam/.tfam, tab-delimited .bim/.map):"
echo "  $GT_PREFIX.bed / .bim / .fam / .tped / .tfam / .map"
echo "Variants:    $(wc -l < "$GT_PREFIX".bim)"
echo "Samples:     $(wc -l < "$GT_PREFIX".fam)"
