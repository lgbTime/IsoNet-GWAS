#!/usr/bin/env bash
# =============================================================================
# run_layer2.sh — Step 5: Layer-2 expression-network expansion
# -----------------------------------------------------------------------------
# IsoNet-GWAS Layer 2 couples the Layer-1 co-function modules with isoform-level
# expression data:
#   * prefilter_expression_matrix.py (optional)
#       - removes low-expression isoforms from the TPM matrix
#   * ISONET_LAYER2.py
#       - expands each AI-confirmed HIGH-confidence Layer-1 module by
#         co-expression (spearman rho >= $CORR_THRESHOLD, FDR <= $FDR_THRESHOLD)
#       - ranks regulators by a combined GWAS + expression score and reports
#         "rising hubs" (isoforms whose rank rises after expansion)
#   * plot_layer2_network.R — renders each expanded module network
#
# Inputs : Layer-1 outputs in $L1_DIR (sig_*_cofunction_{modules,nodes,edges}.tsv,
#          sig_*_network_ai_analysis.json), $EXPRESSION_MATRIX, $GTF, $ANNOTATIONS
# Outputs: $L2_DIR/  (module_*_coexpression_edges.tsv,
#                     module_*_regulatory_ranking.tsv, *_rising_hubs.tsv, plots)
#
# Usage   : source ../config.env && bash run_layer2.sh
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
source "$REPO_ROOT/config.env"

L1_DIR="${L1_DIR:-./layer1_cofunction}"
L2_DIR="${L2_DIR:-./layer2_expression}"
mkdir -p "$L1_DIR" "$L2_DIR"

# Resolve actual Layer-1 output filenames (they embed the phenotype name)
export QQMAN_DIR="${QQMAN_DIR:-$DATA_DIR/qqman}"
MODULES_TSV="${MODULES_TSV:-$(ls -t "$L1_DIR"/sig_*_cofunction_modules.tsv 2>/dev/null | head -1)}"
NODES_TSV="${NODES_TSV:-$(ls -t "$L1_DIR"/sig_*_cofunction_nodes.tsv 2>/dev/null | head -1)}"
EDGES_TSV="${EDGES_TSV:-$(ls -t "$L1_DIR"/sig_*_cofunction_edges.tsv 2>/dev/null | head -1)}"
AI_JSON="${AI_JSON:-$(ls -t "$L1_DIR"/sig_*_network_ai_analysis.json 2>/dev/null | head -1)}"
QQMAN_FILE="${QQMAN_FILE:-$(ls "$QQMAN_DIR"/$QQMAN_PATTERN 2>/dev/null | head -1)}"

for f in "$MODULES_TSV" "$NODES_TSV" "$EDGES_TSV" "$AI_JSON" "$QQMAN_FILE"; do
    if [ -z "$f" ] || [ ! -f "$f" ]; then
        echo "ERROR: missing Layer-1/qqman input: $f" >&2
        echo "       Run run_layer1.sh first (AI analysis required for module confidence)." >&2
        exit 1
    fi
done
echo "  Using Layer-1 inputs:"
echo "    modules: $(basename "$MODULES_TSV")"
echo "    nodes:   $(basename "$NODES_TSV")"
echo "    edges:   $(basename "$EDGES_TSV")"
echo "    AI json: $(basename "$AI_JSON")"
echo "    qqman:   $(basename "$QQMAN_FILE")"

# ---------------------------------------------------------------------------
# Optional: pre-filter the expression matrix to remove low-expression noise
# ---------------------------------------------------------------------------
if [ -n "${PREFILTER_EXPRESSION:-}" ] && [ ! -f "$EXPRESSION_MATRIX.filtered" ]; then
    echo "=== Step 5a: pre-filtering expression matrix ==="
    "$PYTHON" "$SCRIPT_DIR/prefilter_expression_matrix.py" \
        --input "$EXPRESSION_MATRIX" \
        --output "$EXPRESSION_MATRIX.filtered"
    EXPRESSION_MATRIX="$EXPRESSION_MATRIX.filtered"
fi

# ---------------------------------------------------------------------------
# Extract all "high" confidence module IDs from the AI analysis JSON
# ---------------------------------------------------------------------------
HIGH_MODULES=$("$PYTHON" - "$AI_JSON" <<'EOF'
import json, sys
with open(sys.argv[1]) as f:
    data = json.load(f)
high_ids = [m['module_id'] for m in data.get('modules', []) if m.get('confidence') == 'high']
print(' '.join(high_ids))
EOF
)
echo "  High-confidence modules: ${HIGH_MODULES:-<none>}"

if [ -z "$HIGH_MODULES" ]; then
    echo "  ERROR: no high-confidence modules found in $(basename "$AI_JSON")" >&2
    echo "         Rerun Layer-1 with AI enabled, or set MODULES manually." >&2
    exit 1
fi

# Clean any stale module files from previous runs
rm -f "$L2_DIR"/module_*

# ---------------------------------------------------------------------------
# Layer 2: expression expansion of the high-confidence modules
# ---------------------------------------------------------------------------
echo ""
echo "=== Layer 2: Expression expansion (HIGH-confidence modules only) ==="
"$PYTHON" "$SCRIPT_DIR/ISONET_LAYER2.py" \
    --modules-tsv "$MODULES_TSV" \
    --nodes-tsv "$NODES_TSV" \
    --edges-tsv "$EDGES_TSV" \
    --ai-analysis-json "$AI_JSON" \
    --annotations "$ANNOTATIONS" \
    --expression "$EXPRESSION_MATRIX" \
    --gtf "$GTF" \
    --qqman "$QQMAN_FILE" \
    --go-tab "$GO_TAB" \
    --out-prefix "$L2_DIR/out" \
    --ai-confidence high \
    --modules $HIGH_MODULES \
    --min-shared-terms "${MIN_SHARED_TERMS:-2}" \
    --max-expansion "${MAX_EXPANSION:-500}" \
    --expansion-hops "${EXPANSION_HOPS:-1}" \
    --corr-threshold "${CORR_THRESHOLD:-0.65}" \
    --fdr-threshold "${FDR_THRESHOLD:-0.05}" \
    --min-tpm "${MIN_TPM:-0.1}" \
    --gwas-weight "${GWAS_WEIGHT:-0.3}" \
    --expr-weight "${EXPR_WEIGHT:-0.7}" \
    --min-rank-rise "${MIN_RANK_RISE:-10}" \
    --min-gwas-pvalue "${MIN_GWAS_PVALUE:-0.01}" \
    --tf-min-corr "${TF_MIN_CORR:-0.3}" \
    --tf-top-n "${TF_TOP_N:-300}"

# Flatten: strip out_ prefix
for f in "$L2_DIR"/out_module_*; do
    [ -f "$f" ] || continue
    base=$(basename "$f")
    mv "$f" "$L2_DIR/${base#out_}"
done
mv "$L2_DIR"/out_rising_hubs.tsv "$L2_DIR"/rising_hubs.tsv 2>/dev/null || true
rm -f "$L2_DIR"/out_module_*

# ---------------------------------------------------------------------------
# Plot each expanded module / TF sub-module (figures need R; tables do not)
# ---------------------------------------------------------------------------
HAS_R=0
command -v "$RSCRIPT" >/dev/null 2>&1 && HAS_R=1
[ "$HAS_R" -eq 1 ] || echo "  NOTE: Rscript not found — skipping figures (ranking/edge tables still produced)."

if [ "$HAS_R" -eq 1 ]; then
echo ""
echo "=== Plotting expanded modules ==="
for EDGES_FILE in "$L2_DIR"/module_*_coexpression_edges.tsv; do
    [ -f "$EDGES_FILE" ] || continue
    fname=$(basename "$EDGES_FILE")
    mid="${fname#module_}"
    mid="${mid%_coexpression_edges.tsv}"

    M_RANKING="$L2_DIR/module_${mid}_regulatory_ranking.tsv"
    if [ ! -f "$M_RANKING" ]; then
        echo "  WARNING: Missing ranking for $mid"
        continue
    fi

    echo "  -> Plotting $mid..."
    HIGHLIGHT_GENES="${HIGHLIGHT_GENES:-}" LABEL_DENSITY="${LABEL_DENSITY:-dense}" \
    "$RSCRIPT" "$SCRIPT_DIR/plot_layer2_network.R" \
        "$EDGES_FILE" \
        "$M_RANKING" \
        "$NODES_TSV" \
        "$L2_DIR/${mid}" \
        "$mid" \
        "$AI_JSON"
done
fi

echo ""
echo "  $L2_DIR/:"
for EDGES_FILE in "$L2_DIR"/module_*_coexpression_edges.tsv; do
    [ -f "$EDGES_FILE" ] || continue
    fname=$(basename "$EDGES_FILE")
    mid="${fname#module_}"
    mid="${mid%_coexpression_edges.tsv}"
    echo "    ${mid}_coexpression_{network,focused}.{pdf,png,tiff}"
    echo "    ${mid}_evidence_scatter.pdf  ${mid}_kME_barchart.pdf"
done
echo "    rising_hubs.tsv"
echo "  DONE"
