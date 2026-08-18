#!/usr/bin/env bash
# =============================================================================
# 00_master_pipeline.sh — end-to-end IsoNet-GWAS pipeline
# -----------------------------------------------------------------------------
# IsoNet-GWAS: Isoform-PAV GWAS coupled with a two-layer network
#
#   Step 1  (01_isoform_genetic_map)  : isoform PAV -> unique-position
#                                       genetic map -> PLINK dataset
#   Step 2  (02_genotype_plus_marker) : merge extra marker VCFs (gene-body SV,
#                                       GLK PAV, ...) into the PAV genotype
#   Step 3  (03_gwas)                 : EMMAX association per trait
#                                       (phenotype qq-norm -> kinship -> GWAS
#                                        -> .qqman -> QQ/Manhattan plots)
#   Step 4  (04_layer1_cofunction)    : LAYER 1 — significant isoforms ->
#                                       co-function network (GO/PFAM/KEGG/COG)
#                                       + modules + optional AI module calling
#   Step 5  (05_layer2_expression)    : LAYER 2 — expression-based expansion
#                                       of high-confidence modules,
#                                       regulatory ranking, rising hubs
#
# Usage:
#   ./00_master_pipeline.sh            # run steps 1-5
#   ./00_master_pipeline.sh 3 5        # run only steps 3,4,5
#
# Before running: edit config.env (tool paths, data locations, AI keys).
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"
source ./config.env

FROM_STEP=1
TO_STEP=5
for arg in "$@"; do
    case "$arg" in
        [1-5]) if [ "$FROM_STEP" -eq 1 ] && [ "$TO_STEP" -eq 5 ]; then
                   FROM_STEP="$arg"; TO_STEP="$arg";
               else
                   TO_STEP="$arg";
               fi ;;
        *) echo "Unknown argument: $arg (expected step numbers 1-5)" >&2; exit 1 ;;
    esac
done
if [ "$FROM_STEP" -lt 1 ] || [ "$TO_STEP" -gt 5 ] || [ "$FROM_STEP" -gt "$TO_STEP" ]; then
    echo "ERROR: invalid step range ($FROM_STEP..$TO_STEP); steps must be within 1..5" >&2
    exit 1
fi

echo "================================================================"
echo " IsoNet-GWAS master pipeline  (steps $FROM_STEP-$TO_STEP)"
echo "================================================================"

# ---------------------------------------------------------------------------
# Dependency checks
# ---------------------------------------------------------------------------
need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not found in PATH (see config.env)" >&2; exit 1; }; }
need "$PYTHON"
[ "$FROM_STEP" -le 2 ] && [ "$TO_STEP" -ge 1 ] && need "$PLINK"
[ "$FROM_STEP" -le 3 ] && [ "$TO_STEP" -ge 3 ] && need "$EMMAX" && need "$EMMAX_KIN"
if [ "$TO_STEP" -ge 4 ] && ! command -v "$RSCRIPT" >/dev/null 2>&1; then
    echo "ERROR: Rscript not found in PATH — R is required for the pipeline" >&2
    echo "       (phenotype qq-norm, QQ/Manhattan and network figures)." >&2
    echo "       Set RSCRIPT in config.env to your Rscript binary." >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Step 1 — isoform genetic map + PLINK dataset
# ---------------------------------------------------------------------------
run_step1() {
    echo ""
    echo "########## STEP 1: Isoform genetic map -> PLINK dataset ##########"
    bash "$SCRIPT_DIR/01_isoform_genetic_map/bed2map/get_uniq_pos.sh"
    bash "$SCRIPT_DIR/01_isoform_genetic_map/pav2ped/plink.sh"
}

# ---------------------------------------------------------------------------
# Step 2 — plus-marker genotype
# ---------------------------------------------------------------------------
run_step2() {
    echo ""
    echo "########## STEP 2: Merge extra marker VCFs into genotype ##########"
    if [ ! -f "$GENE_BODY_VCF" ] && [ ! -f "$GLK_PAV_VCF" ]; then
        echo "NOTE: neither GENE_BODY_VCF nor GLK_PAV_VCF exists — skipping Step 2."
        echo "      (The Step-1 PAV genotype ($PAV_PREFIX) is used for GWAS.)"
        export GT_PREFIX="$PAV_PREFIX"
        return 0
    fi
    bash "$SCRIPT_DIR/02_genotype_plus_marker/run.sh"
}

# ---------------------------------------------------------------------------
# Step 3 — GWAS
# ---------------------------------------------------------------------------
run_step3() {
    echo ""
    echo "########## STEP 3: EMMAX GWAS per trait ##########"
    if [ ! -f "$DATA_DIR/map/genetic_map.map" ]; then
        echo "ERROR: genetic map not found: $DATA_DIR/map/genetic_map.map" >&2
        echo "       Run Step 1 first." >&2
        exit 1
    fi
    bash "$SCRIPT_DIR/03_gwas/gwas.sh"
}

# ---------------------------------------------------------------------------
# Step 4 — Layer 1 co-function network
# ---------------------------------------------------------------------------
run_step4() {
    echo ""
    echo "########## STEP 4: Layer-1 co-function network ##########"
    # Point the network step at the GWAS output directory produced in Step 3
    export QQMAN_DIR="${QQMAN_DIR:-$SCRIPT_DIR/03_gwas/gwasout}"
    bash "$SCRIPT_DIR/04_layer1_cofunction_network/run_layer1.sh"
}

# ---------------------------------------------------------------------------
# Step 5 — Layer 2 expression network
# ---------------------------------------------------------------------------
run_step5() {
    echo ""
    echo "########## STEP 5: Layer-2 expression network ##########"
    bash "$SCRIPT_DIR/05_layer2_expression_network/run_layer2.sh"
}

# ---------------------------------------------------------------------------
# Execute selected steps
# ---------------------------------------------------------------------------
[ "$FROM_STEP" -le 1 ] && [ "$TO_STEP" -ge 1 ] && run_step1
[ "$FROM_STEP" -le 2 ] && [ "$TO_STEP" -ge 2 ] && run_step2
[ "$FROM_STEP" -le 3 ] && [ "$TO_STEP" -ge 3 ] && run_step3
[ "$FROM_STEP" -le 4 ] && [ "$TO_STEP" -ge 4 ] && run_step4
[ "$FROM_STEP" -le 5 ] && [ "$TO_STEP" -ge 5 ] && run_step5

echo ""
echo "================================================================"
echo " IsoNet-GWAS pipeline finished (steps $FROM_STEP-$TO_STEP)."
echo "================================================================"
