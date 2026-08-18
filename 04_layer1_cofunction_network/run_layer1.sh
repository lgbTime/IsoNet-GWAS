#!/usr/bin/env bash
# =============================================================================
# run_layer1.sh — Step 4: Layer-1 co-function network
# -----------------------------------------------------------------------------
# IsoNet-GWAS Layer 1 couples the GWAS results with a co-function network:
#   * ISONET_LAYER1.py
#       - collects significant isoforms from all *.pheno.qqman files
#       - maps each isoform to its parent gene (StringTie GTF / reference GFF)
#       - annotates via eggNOG-mapper output (GO/PFAM/KEGG/COG/description)
#       - builds a co-function network (edges = shared functional terms)
#       - calls network modules (community detection)
#       - optionally uses an LLM (DeepSeek/OpenAI) to name and prioritize
#         regulatory modules (set DEEPSEEK_API_KEY in config.env)
#   * build_v1_concentric.py + plot_layer1_network.R
#       - renders the V1 concentric-circle network figure
#
# Inputs : $QQMAN_DIR/*.pheno.qqman   (GWAS results from Step 3)
#          $GTF, $GFF, $ANNOTATIONS, $GO_TAB
# Outputs: $L1_DIR/  (sig_*_cofunction_{nodes,edges,modules}.tsv,
#                     *_network_*.pdf/png, *_network_ai_analysis.json, ...)
#
# Usage   : source ../config.env && bash run_layer1.sh
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
source "$REPO_ROOT/config.env"

L1_DIR="${L1_DIR:-./layer1_cofunction}"
export QQMAN_DIR="${QQMAN_DIR:-$DATA_DIR/qqman}"
mkdir -p "$L1_DIR"

if ! ls "$QQMAN_DIR"/$QQMAN_PATTERN >/dev/null 2>&1; then
    echo "ERROR: no files matching '$QQMAN_PATTERN' in QQMAN_DIR=$QQMAN_DIR" >&2
    exit 1
fi

# Assemble optional AI arguments (skip entirely when no API key is set)
AI_ARGS=()
if [ -n "${DEEPSEEK_API_KEY:-}" ]; then
    AI_ARGS=(--api-key "$DEEPSEEK_API_KEY" \
             --model "$DEEPSEEK_MODEL" \
             --api-base "$DEEPSEEK_API_BASE" \
             --thinking)
    echo "  AI enrichment enabled (model: $DEEPSEEK_MODEL)"
else
    echo "  NOTE: DEEPSEEK_API_KEY not set — skipping AI module analysis."
    echo "        Network construction still runs."
fi

# R is only needed for figures — network tables/TSVs are produced without it
HAS_R=0
command -v "$RSCRIPT" >/dev/null 2>&1 && HAS_R=1
[ "$HAS_R" -eq 1 ] || echo "  NOTE: Rscript not found — skipping figures (network tables still produced)."

# ═══════════════ Layer 1: GWAS -> co-function network + AI ═══════════════
echo "=== Layer 1: GWAS + co-function + AI ==="
PLOT_ARGS=()
[ "$HAS_R" -eq 1 ] && PLOT_ARGS=(--plot-network)
"$PYTHON" "$SCRIPT_DIR/ISONET_LAYER1.py" \
    --qqman-dir "$QQMAN_DIR" \
    --qqman-pattern "$QQMAN_PATTERN" \
    --gtf "$GTF" \
    --gff "$GFF" \
    --annotations "$ANNOTATIONS" \
    --go-tab "$GO_TAB" \
    --pvalue-threshold "$P_VALUE_THRESHOLD" \
    --out-prefix "$L1_DIR/sig" \
    --build-network \
    --network-top-n 300 \
    --network-min-shared 2 \
    --network-min-module-size 3 \
    "${PLOT_ARGS[@]}" \
    "${AI_ARGS[@]}"

# ═══════════════ V1 concentric visualization ═══════════════
# (full + cofunction + cofunction_orphans views only)
QQMAN_FILE="${QQMAN_FILE:-$(ls "$QQMAN_DIR"/$QQMAN_PATTERN | head -1)}"
echo "  -> V1 concentric (based on $(basename "$QQMAN_FILE"))..."
"$PYTHON" "$SCRIPT_DIR/build_v1_concentric.py" \
    "$QQMAN_FILE" "$ANNOTATIONS" "$GFF" "$GO_TAB" "$L1_DIR"

if [ "$HAS_R" -eq 1 ]; then
    PLOT_VIEWS="${PLOT_VIEWS:-full,cofunction,cofunction_orphans}" \
    "$RSCRIPT" "$SCRIPT_DIR/plot_layer1_network.R" \
        "$L1_DIR/v1_nodes.tsv" \
        "$L1_DIR/v1_edges.tsv" \
        "$L1_DIR/sig_concentric"
else
    echo "  (skipped: v1 concentric figure needs Rscript)"
fi

echo ""
echo "  $L1_DIR/:"
echo "    sig_*_cofunction_{nodes,edges,modules}.tsv"
echo "    sig_*_network_cofunction{,_focused}.{pdf,png}"
echo "    sig_*_network_ai_analysis.json  (when AI enabled)"
echo "    sig_concentric_network_{full,cofunction,cofunction_orphans}.{pdf,png}"
echo "  DONE"
