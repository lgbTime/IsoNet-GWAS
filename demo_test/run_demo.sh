#!/usr/bin/env bash
# =============================================================================
# run_demo.sh — GLK demo: gwasor GWAS → Layer-1 co-function → Layer-2 expression
# -----------------------------------------------------------------------------
# Runs the complete IsoNet-GWAS flow on the small GLK demo dataset (134 lettuce
# accessions, isoform-PAV genotype + the GLK 1-bp deletion marker) without
# needing any API key:
#
#   1. gwasor pri_gwas   : EMMAX association for the chlorophyll a+b trait
#      (phenotype qq-norm → kinship → association → .ps → .qqman → plots)
#   2. ISONET_LAYER1.py : significant isoforms → co-function network
#      (GO/PFAM/KEGG/COG edges) → network modules (+ V1 concentric figure)
#   3. ISONET_LAYER2.py : expression-based expansion of the modules
#      (co-expression + regulatory ranking + rising hubs)
#
# Optional: set DEEPSEEK_API_KEY to also run the AI module-calling step.
# Set DEMO_MODULES="M1 M3" to restrict Layer-2 expansion to specific modules
# (default: all modules found by Layer-1).
#
# Usage   : bash run_demo.sh
# Outputs : output/gwasout_all_chlorophy_a+b/   — GWAS results + QQ/Manhattan
#           output/layer1/                      — co-function network + modules
#           output/layer2/                      — expression expansion + plots
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$SCRIPT_DIR"
source "$REPO_ROOT/config.env"

DEMO_DATA="$SCRIPT_DIR/data"
OUT="$SCRIPT_DIR/output"
mkdir -p "$OUT/layer1" "$OUT/layer2"

# Sanity check: demo data present?
[ -f "$DEMO_DATA/addmode_Only4GLK.bed" ] || {
    echo "ERROR: demo data missing — run: bash build_demo_data.sh" >&2; exit 1; }

# ---------------------------------------------------------------------------
# Step 1 — GWAS (gwasor pri_gwas, EMMAX)
# ---------------------------------------------------------------------------
echo ""
echo "═══════════════ DEMO STEP 1/3: GWAS (gwasor pri_gwas) ═══════════════"
"$PYTHON" "$REPO_ROOT/03_gwas/gwasor.py" pri_gwas \
    -b   "$DEMO_DATA/addmode_Only4GLK" \
    -p   "$DEMO_DATA/all_chlorophy_a+b" \
    -t   no \
    -sig "$SIG" \
    -o   "$OUT/gwasout_all_chlorophy_a+b" \
    -map "$DEMO_DATA/genetic_map.map"

# The network steps expect the *.pheno.qqman naming (phenotype-name prefix)
cp "$OUT/gwasout_all_chlorophy_a+b/all_chlorophy_a+b.qqman" "$OUT/chlorophyll.pheno.qqman"

# ---------------------------------------------------------------------------
# Step 2 — Layer 1: co-function network
# ---------------------------------------------------------------------------
echo ""
echo "═══════════════ DEMO STEP 2/3: Layer-1 co-function network ═══════════════"
HAS_R=0; command -v "$RSCRIPT" >/dev/null 2>&1 && HAS_R=1
PLOT_ARGS=(); [ "$HAS_R" -eq 1 ] && PLOT_ARGS=(--plot-network)
AI_ARGS=()
if [ -n "${DEEPSEEK_API_KEY:-}" ]; then
    AI_ARGS=(--api-key "$DEEPSEEK_API_KEY" --model "$DEEPSEEK_MODEL" \
             --api-base "$DEEPSEEK_API_BASE" --thinking)
fi

"$PYTHON" "$REPO_ROOT/04_layer1_cofunction_network/ISONET_LAYER1.py" \
    --qqman-dir "$OUT" \
    --qqman-pattern "chlorophyll.pheno.qqman" \
    --gtf "$DEMO_DATA/pop240GTF2ref.combined.gtf" \
    --gff "$DEMO_DATA/Lsat_Salinas_v11_genomic.gff" \
    --annotations "$DEMO_DATA/isoform_and_gene.emapper.annotations" \
    --go-tab "$DEMO_DATA/go.tab" \
    --pvalue-threshold "$P_VALUE_THRESHOLD" \
    --out-prefix "$OUT/layer1/sig" \
    --build-network \
    --network-top-n 300 \
    --network-min-shared 2 \
    --network-min-module-size 3 \
    "${PLOT_ARGS[@]}" \
    "${AI_ARGS[@]}"

# ---------------------------------------------------------------------------
# Step 3 — Layer 2: expression expansion (no AI json needed: explicit modules)
# ---------------------------------------------------------------------------
echo ""
echo "═══════════════ DEMO STEP 3/3: Layer-2 expression network ═══════════════"
MODULES_TSV="$(ls -t "$OUT"/layer1/sig_*_cofunction_modules.tsv | head -1)"
NODES_TSV="$(ls -t "$OUT"/layer1/sig_*_cofunction_nodes.tsv | head -1)"
EDGES_TSV="$(ls -t "$OUT"/layer1/sig_*_cofunction_edges.tsv | head -1)"
ALL_MODULES=$(awk 'NR>1{print $1}' "$MODULES_TSV" | tr '\n' ' ')
DEMO_MODULES="${DEMO_MODULES:-$ALL_MODULES}"
echo "  Layer-1 modules: $ALL_MODULES"
echo "  Expanding:       $DEMO_MODULES"

"$PYTHON" "$REPO_ROOT/05_layer2_expression_network/ISONET_LAYER2.py" \
    --modules-tsv "$MODULES_TSV" \
    --nodes-tsv "$NODES_TSV" \
    --edges-tsv "$EDGES_TSV" \
    --annotations "$DEMO_DATA/isoform_and_gene.emapper.annotations" \
    --expression "$DEMO_DATA/Pop_isoform_TPM_by_compareGTF.txt" \
    --gtf "$DEMO_DATA/pop240GTF2ref.combined.gtf" \
    --qqman "$OUT/chlorophyll.pheno.qqman" \
    --go-tab "$DEMO_DATA/go.tab" \
    --out-prefix "$OUT/layer2/out" \
    --modules $DEMO_MODULES \
    --min-shared-terms 2 \
    --max-expansion 500 \
    --expansion-hops 1 \
    --corr-threshold 0.65 \
    --fdr-threshold 0.05 \
    --min-tpm 0.1 \
    --gwas-weight 0.3 \
    --expr-weight 0.7 \
    --min-rank-rise 10 \
    --min-gwas-pvalue 0.01 \
    --tf-min-corr 0.3 \
    --tf-top-n 300

# Flatten: strip out_ prefix (keeps run_layer2.sh naming convention)
for f in "$OUT"/layer2/out_module_*; do
    [ -f "$f" ] || continue
    mv "$f" "$OUT/layer2/${f##*/out_}"
done
mv "$OUT"/layer2/out_rising_hubs.tsv "$OUT"/layer2/rising_hubs.tsv 2>/dev/null || true

# ---------------------------------------------------------------------------
# Plot each expanded module (figures need R; tables do not)
# ---------------------------------------------------------------------------
if [ "$HAS_R" -eq 1 ]; then
    echo ""
    echo "=== Plotting expanded modules ==="
    for EDGES_FILE in "$OUT"/layer2/module_*_coexpression_edges.tsv; do
        [ -f "$EDGES_FILE" ] || continue
        fname=$(basename "$EDGES_FILE")
        mid="${fname#module_}"
        mid="${mid%_coexpression_edges.tsv}"
        M_RANKING="$OUT/layer2/module_${mid}_regulatory_ranking.tsv"
        [ -f "$M_RANKING" ] || continue
        echo "  -> Plotting $mid..."
        HIGHLIGHT_GENES="${HIGHLIGHT_GENES:-}" LABEL_DENSITY="${LABEL_DENSITY:-dense}" \
        "$RSCRIPT" "$REPO_ROOT/05_layer2_expression_network/plot_layer2_network.R" \
            "$EDGES_FILE" "$M_RANKING" "$NODES_TSV" \
            "$OUT/layer2/${mid}" "$mid" "NONE"
    done
fi

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
echo ""
echo "════════════════════════ DEMO DONE ════════════════════════"
echo "  GWAS:    $OUT/gwasout_all_chlorophy_a+b/  (.ps/.qqman/.man.jpeg/.qq.jpeg)"
echo "  Layer 1: $OUT/layer1/sig_chlorophyll_cofunction_{nodes,edges,modules}.tsv"
echo "  Layer 2: $OUT/layer2/ (module_*_regulatory_ranking.tsv, rising_hubs.tsv)"
if [ "$HAS_R" -eq 1 ]; then
    echo "  Figures: $OUT/layer1/*_network_*.pdf   $OUT/layer2/module_*_coexpression_network.pdf"
else
    echo "  (Rscript not found — figures skipped, tables still produced)"
fi
