#!/usr/bin/env python3
"""
ISONET_LAYER2.py : Expression-based expansion and re-ranking
                                       of co-function modules.

For each co-function module from ISONET_LAYER1:
  1. Extract characterizing functional terms (GO / PFAM / KEGG / COG / Description)
  2. Expand module by searching the FULL eggNOG annotation space for isoforms
     sharing those terms : capturing functionally related genes missed by GWAS
  3. Build a co-expression network from TPM data within the expanded module
  4. Compute multiple centrality metrics (degree, betweenness, kME)
  5. Integrate GWAS significance + expression centrality into a combined score
  6. Identify "rising hubs" : genes that rank low in GWAS but high in expression
     centrality, flagging them as candidate regulators

Rationale:
    GWAS association strength ? causal regulatory importance.  A gene may be a
    weak GWAS hit (or miss the significance threshold entirely) yet be the
    central regulatory hub in the co-expression network of its functional
    module.  Layer-2 expression analysis recovers these hidden regulators by
    triangulating functional annotation + co-expression + GWAS signal.

Output files (per module):
    {prefix}_module_{name}_expanded_nodes.tsv         : expanded gene set
    {prefix}_module_{name}_coexpression_edges.tsv     : significant co-expression edges
    {prefix}_module_{name}_regulatory_ranking.tsv     : combined-score ranking
    {prefix}_rising_hubs.tsv                          : cross-module rising hub report
    {prefix}_module_{name}_coexpression_network.pdf   : (optional) R plot

Usage:
    python3 ISONET_LAYER2.py \\
        --modules-tsv out_GLK_merge_qqman/chrolophyll_significant_chlorophyll_cofunction_modules.tsv \\
        --nodes-tsv   out_GLK_merge_qqman/chrolophyll_significant_chlorophyll_cofunction_nodes.tsv \\
        --edges-tsv   out_GLK_merge_qqman/chrolophyll_significant_chlorophyll_cofunction_edges.tsv \\
        --annotations isoform_and_gene.emapper.annotations \\
        --expression  Pop_isoform_TPM_by_compareGTF.txt \\
        --gtf         pop240GTF2ref.combined.gtf \\
        --qqman       chlorophyll.pheno.qqman \\
        --out-prefix  layer2_chlorophyll
"""

import argparse
import json
import os
import re
import sys
import time
import urllib.parse
from collections import Counter, defaultdict
from itertools import combinations

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import spearmanr
import networkx as nx


# ---------------------------------------------------------------------------
# Generic GO / KEGG term filters (same as ISONET_LAYER1.py)
# ---------------------------------------------------------------------------

_GENERIC_GO_TERMS = {
    'catalytic activity', 'binding', 'atp binding', 'protein binding',
    'dna binding', 'metal ion binding', 'transferase activity',
    'hydrolase activity', 'nucleic acid binding', 'rna binding',
    'nucleotide binding', 'transferase activity, transferring '
        'phosphorus-containing groups',
    'kinase activity', 'oxidoreductase activity', 'lyase activity',
    'isomerase activity', 'ligase activity', 'transporter activity',
    'molecular_function', 'biological_process', 'cellular_component',
    'intracellular anatomical structure', 'obsolete cell', 'cytoplasm',
    'nucleus', 'plasma membrane', 'membrane', 'chloroplast', 'plastid',
    'mitochondrion', 'cytosol', 'endoplasmic reticulum',
    'golgi apparatus', 'extracellular region', 'cell wall',
    'obsolete intracellular', 'intracellular', 'organelle',
    'protein-containing complex', 'cellular anatomical entity',
    'metabolic process', 'biosynthetic process', 'catabolic process',
    'cellular process', 'biological regulation', 'regulation of '
        'biological process',
    'response to stress', 'transport', 'signal transduction',
    'protein phosphorylation', 'phosphorylation',
    'oxidation-reduction process', 'proteolysis',
    'defense response', 'response to stimulus',
    'nucleobase-containing compound metabolic process',
    'dna metabolic process', 'rna metabolic process',
    'dna repair', 'dna replication', 'transcription, dna-templated',
    'translation', 'protein folding', 'protein transport',
    'intracellular protein transport', 'vesicle-mediated transport',
    'lipid metabolic process', 'carbohydrate metabolic process',
    'amino acid metabolic process', 'protein metabolic process',
    'response to acid chemical', 'response to salt stress',
    'function unknown', 'general function prediction only',
    'reproduction', 'developmental process involved in reproduction',
    'multicellular organismal process', 'growth',
    'anatomical structure development', 'cell differentiation',
    'post-embryonic development', 'regulation of transcription, dna-templated',
    'regulation of dna-templated transcription',
}

_GENERIC_COG_EDGE_TERMS = {
    '[K] Transcription',
    '[T] Signal transduction mechanisms',
    '[R] General function prediction only',
    '[S] Function unknown',
}

_BROAD_KEGG_IDS = {
    'map01100', 'ko01100',     # Metabolic pathways (global overview)
    'map01110', 'ko01110',     # Biosynthesis of secondary metabolites (broad)
}

# Vague description patterns (too generic for co-function signal)
_VAGUE_DESC_RE = re.compile(
    r'^(uncharacteri[sz]ed|hypothetical|predicted|probable|putative|'
    r'expressed|unknown|unnamed|orf|cds|partial|'
    r'encoded by|belongs to the|contain|domain|family|'
    r'protein of unknown|duf\d+)',
    re.IGNORECASE)

# ---------------------------------------------------------------------------
# Transcription factor module detection
# ---------------------------------------------------------------------------

# PFAM domains that indicate a transcription factor module
_TF_PFAM_TERMS = {
    'Myb_DNA-binding', 'Myb_CC_LHEQLE', 'Myb_CC', 'Myb',
    'bZIP_1', 'bZIP_2', 'bZIP_Maf', 'bZIP',
    'WRKY', 'AP2', 'ERF', 'GATA', 'MADS', 'TCP', 'NAC', 'GRAS',
    'Homeobox', 'Homeodomain', 'HLH', 'bHLH', 'bHLH-MYC',
    'Zinc finger', 'zf-C2H2', 'zf-C3HC4', 'zf-Dof',
    'E2F_TDP', 'CBFB_NFYA', 'NF-YB', 'NF-YC',
    'TF_bZIP', 'TF_Otx', 'TF_AP-2',
    'HMG_box', 'SRF-TF', 'Ets', 'Forkhead', 'PAX',
    'SBP', 'ARID', 'CSD', 'CBF', 'CPP', 'DBB', 'EIL',
    'FAR1', 'G2-like', 'GeBP', 'HRT', 'HSF_DNA-bind',
    'LBD', 'LFY', 'LSD', 'M-type_MADS', 'MIKC_MADS',
    'NF-X1', 'NF-YA', 'NF-YB', 'NF-YC', 'NZF',
    'RWP-RK', 'S1FA', 'SAP', 'STAT', 'T-box', 'TIG',
    'Trihelix', 'VOZ', 'Whirly', 'YABBY', 'ZF-HD',
    'B3', 'ARF', 'Auxin_resp', 'AUX_IAA',
}


def _is_tf_module(module_terms):
    """Return True if the module is dominated by transcription factor terms."""
    pfam_terms = {t for t, _ in module_terms.get('pfam', [])}
    tf_hits = pfam_terms & _TF_PFAM_TERMS
    return len(tf_hits) >= 1


def _load_gene_expression_matrix(expression_path, expr_available_ids=None):
    """Load and aggregate the full expression matrix to gene level.

    Reads the isoform-level TPM matrix and sums isoforms into gene-level
    TPMs.  Caches the result on the function for reuse across modules.

    Returns
    -------
    gene_expr : pd.DataFrame
        Rows = gene IDs (LOC*), columns = sample names.  TPMs summed across
        all measured isoforms of each gene.
    """
    # Cache on the function attribute (singleton pattern)
    cache_attr = '_gene_expr_cache'
    if hasattr(_load_gene_expression_matrix, cache_attr):
        return getattr(_load_gene_expression_matrix, cache_attr)

    print(f"  [TF-TARGET] Loading full gene-level expression matrix …")
    t0 = time.time()

    gene_tpm = defaultdict(lambda: np.zeros(0))
    sample_names = None

    with open(expression_path, 'r') as fh:
        header = fh.readline().rstrip('\n').split('\t')
        sample_names = header[1:]
        n_samples = len(sample_names)

        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 2:
                continue
            expr_key = parts[0]
            if ':' in expr_key:
                loc_id = expr_key.split(':', 1)[0]
            else:
                loc_id = expr_key

            # Skip genes not in the pre-scanned expression IDs (optional filter)
            if expr_available_ids and loc_id not in expr_available_ids:
                continue

            try:
                values = np.array([float(v) for v in parts[1:]], dtype=np.float32)
            except ValueError:
                values = np.zeros(n_samples, dtype=np.float32)

            if loc_id not in gene_tpm:
                gene_tpm[loc_id] = values
            else:
                gene_tpm[loc_id] += values

    # Build DataFrame
    gene_ids = sorted(gene_tpm.keys())
    data = np.array([gene_tpm[g] for g in gene_ids], dtype=np.float32)
    gene_expr = pd.DataFrame(data, index=gene_ids, columns=sample_names)

    elapsed = time.time() - t0
    print(f"  [TF-TARGET]   {len(gene_ids)} genes × {n_samples} samples "
          f"({elapsed:.1f}s)")

    setattr(_load_gene_expression_matrix, cache_attr, gene_expr)
    return gene_expr


def find_tf_targets(expression_path, tf_seed_genes, expr_available_ids=None,
                    top_n=300, min_corr=0.3):
    """Find candidate target genes whose expression correlates with TF seeds.

    For each TF seed that is present in the expression matrix, compute the
    Pearson correlation of its expression vector against every other gene.
    Returns the union of top-N correlated genes across all TF seeds.

    This captures genes whose transcript levels track with a TF regulator
    in vivo — the functional targets that annotation alone would miss.

    Parameters
    ----------
    expression_path : str
        Path to the TPM expression matrix.
    tf_seed_genes : set of str
        Gene IDs (LOC*) of TF seeds whose targets we want to find.
    expr_available_ids : set or None
        Pre-scanned expression IDs (optional filter).
    top_n : int
        Maximum targets to return per TF seed.
    min_corr : float
        Minimum absolute Pearson r for a gene to be considered a target.

    Returns
    -------
    targets : set of str
        Gene IDs (LOC*) whose expression correlates with at least one TF seed.
    tf_target_detail : dict
        tf_seed → list of (target_gene, pearson_r) sorted by |r| descending.
    """
    gene_expr = _load_gene_expression_matrix(expression_path, expr_available_ids)

    # Map TF seeds to their rows in the expression matrix
    tf_seeds_in_expr = {g for g in tf_seed_genes if g in gene_expr.index}
    if not tf_seeds_in_expr:
        print(f"  [TF-TARGET] WARNING: no TF seeds found in expression matrix")
        return set(), {}

    # Also try to find TF seeds via LOC prefix match (gene may be LOC111886262
    # but expression key may differ)
    if len(tf_seeds_in_expr) == 0:
        for seed in tf_seed_genes:
            if seed.startswith('LOC'):
                matches = [g for g in gene_expr.index if g == seed]
                tf_seeds_in_expr.update(matches)

    if not tf_seeds_in_expr:
        return set(), {}

    all_targets = set()
    tf_target_detail = {}

    # Pre-compute expression matrix as numpy for fast vectorized correlation
    expr_mat = gene_expr.values.T  # samples × genes
    gene_list = list(gene_expr.index)
    gene_to_idx = {g: i for i, g in enumerate(gene_list)}
    n_samples = expr_mat.shape[0]

    print(f"  [TF-TARGET] Correlating {len(tf_seeds_in_expr)} TF seed(s) "
          f"against {len(gene_list)} genes (|r| ≥ {min_corr}) …")
    t0 = time.time()

    for tf_gene in sorted(tf_seeds_in_expr):
        if tf_gene not in gene_to_idx:
            continue

        tf_idx = gene_to_idx[tf_gene]
        tf_vec = expr_mat[:, tf_idx].astype(np.float64)

        # Skip low-variance TF expression (no signal)
        if np.std(tf_vec) < 1e-6:
            print(f"  [TF-TARGET]   {tf_gene}: zero-variance, skipping")
            continue

        # Vectorized Pearson correlation against all genes
        # r = dot(x - mean(x), y - mean(y)) / (n * std(x) * std(y))
        tf_centered = tf_vec - tf_vec.mean()
        tf_norm = np.sqrt(np.dot(tf_centered, tf_centered))

        # Compute all correlations at once using matrix-vector product
        all_centered = expr_mat.astype(np.float64) - expr_mat.mean(axis=0)
        all_norms = np.sqrt((all_centered ** 2).sum(axis=0))
        cov = np.dot(tf_centered, all_centered)
        # Avoid division by zero for constant genes
        denom = tf_norm * all_norms
        denom[denom < 1e-12] = 1.0
        r = cov / denom

        # Find top-N by absolute correlation
        abs_r = np.abs(r)
        # Exclude self-correlation
        abs_r[tf_idx] = 0.0

        top_indices = np.argsort(abs_r)[-top_n:][::-1]
        targets = []
        for idx in top_indices:
            if abs_r[idx] >= min_corr:
                targets.append((gene_list[idx], float(r[idx])))

        tf_target_detail[tf_gene] = targets
        target_genes = {t[0] for t in targets}
        all_targets.update(target_genes)

        print(f"  [TF-TARGET]   {tf_gene}: {len(targets)} targets "
              f"(|r| range: {min(abs_r[top_indices[:len(targets)]]):.2f}–"
              f"{max(abs_r[top_indices[:len(targets)]]) if targets else 0:.2f})")

    elapsed = time.time() - t0
    print(f"  [TF-TARGET]   {len(all_targets)} unique target genes "
          f"across {len(tf_seeds_in_expr)} TF(s) ({elapsed:.1f}s)")

    return all_targets, tf_target_detail


def _is_generic_go(term):
    """Check if a GO term description is in the generic filter set."""
    return term.strip().lower() in _GENERIC_GO_TERMS


def _is_broad_kegg(term):
    """Check if a KEGG pathway term is a broad reference map."""
    if not term.startswith('KEGG:'):
        return False
    rest = term[5:]
    colon = rest.find(':')
    if colon == -1:
        return rest in _BROAD_KEGG_IDS
    return rest[:colon] in _BROAD_KEGG_IDS


def _is_vague_description(desc):
    """Check if an eggNOG Description string is too vague for co-function signal."""
    return bool(_VAGUE_DESC_RE.search(desc))


# ---------------------------------------------------------------------------
# Annotation indexing (term ? query_id mapping)
# ---------------------------------------------------------------------------

def build_annotation_index(annotations_path):
    """Parse full eggNOG-mapper annotations and build term?isoform indices.

    Builds separate indices for PFAM domains, GO term descriptions, COG
    categories, KEGG pathway IDs, and non-vague Descriptions.

    Parameters
    ----------
    annotations_path : str
        Path to eggNOG-mapper annotations file.

    Returns
    -------
    index : dict
        Keys: 'pfam', 'go', 'cog', 'kegg', 'desc' ? each maps term?set of query_ids.
    query_info : dict
        query_id ? {pfam, go, cog, kegg, desc} (raw terms for each query).
    """
    print(f"[INDEX] Building annotation term?gene index from: {annotations_path}")
    t0 = time.time()

    idx = {
        'pfam': defaultdict(set),
        'go': defaultdict(set),
        'cog': defaultdict(set),
        'kegg': defaultdict(set),
        'desc': defaultdict(set),
    }
    query_info = {}

    with open(annotations_path, 'r') as fh:
        # Skip comment lines, find header
        header = None
        for line in fh:
            if line.startswith('##'):
                continue
            header = line.strip().lstrip('#').split('\t')
            break

        if header is None:
            raise ValueError("Could not find header in annotation file")

        # Map column names to indices
        col_map = {name: i for i, name in enumerate(header)}
        qid_col = col_map.get('query', 0)
        desc_col = col_map.get('Description', 7)
        go_col = col_map.get('GOs', 9)
        kegg_col = col_map.get('KEGG_Pathway', 12)
        pfam_col = col_map.get('PFAMs', 20)
        cog_col = col_map.get('COG_category', 6)

        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < max(pfam_col, go_col, kegg_col, cog_col) + 1:
                continue

            qid = parts[qid_col]
            info = {'pfam': [], 'go': [], 'cog': [], 'kegg': [], 'desc': ''}

            # PFAM
            pfam_raw = parts[pfam_col] if pfam_col < len(parts) else ''
            if pfam_raw and pfam_raw != '-':
                for t in pfam_raw.split(','):
                    t = t.strip()
                    if t:
                        idx['pfam'][t].add(qid)
                        info['pfam'].append(t)

            # GO
            go_raw = parts[go_col] if go_col < len(parts) else ''
            if go_raw and go_raw != '-':
                # GO IDs are comma-separated; we store them as-is for later
                # resolution to descriptions via go.tab.  Here we index by GO ID.
                for t in go_raw.split(','):
                    t = t.strip()
                    if t and t != '-':
                        idx['go'][t].add(qid)
                        info['go'].append(t)

            # COG
            cog_raw = parts[cog_col] if cog_col < len(parts) else ''
            if cog_raw and cog_raw != '-':
                for t in cog_raw.split(','):
                    t = t.strip()
                    if t and t != '-':
                        idx['cog'][t].add(qid)
                        info['cog'].append(t)

            # KEGG Pathway
            kegg_raw = parts[kegg_col] if kegg_col < len(parts) else ''
            if kegg_raw and kegg_raw != '-':
                for t in kegg_raw.split(','):
                    t = t.strip()
                    if t and t != '-':
                        idx['kegg'][t].add(qid)
                        info['kegg'].append(t)

            # Description (decoded)
            desc_raw = parts[desc_col] if desc_col < len(parts) else ''
            if desc_raw and desc_raw != '-':
                desc_decoded = urllib.parse.unquote(desc_raw)
                if not _is_vague_description(desc_decoded):
                    # Index by first 120 chars as key
                    desc_key = desc_decoded[:120].strip().lower()
                    if desc_key:
                        idx['desc'][desc_key].add(qid)
                        info['desc'] = desc_decoded

            query_info[qid] = info

    elapsed = time.time() - t0
    n_queries = len(query_info)
    n_pfam = len(idx['pfam'])
    n_go = len(idx['go'])
    n_cog = len(idx['cog'])
    n_kegg = len(idx['kegg'])
    n_desc = len(idx['desc'])
    print(f"[INDEX]   {n_queries:,} query entries indexed in {elapsed:.1f}s")
    print(f"[INDEX]   PFAM: {n_pfam:,} terms  |  GO: {n_go:,} IDs  |  "
          f"COG: {n_cog:,}  |  KEGG: {n_kegg:,}  |  DESC: {n_desc:,} keys")
    return idx, query_info


# ---------------------------------------------------------------------------
# Module characterizing term extraction
# ---------------------------------------------------------------------------

def _process_one_expansion(mid, mod_name, expanded, mod_isos_in_nodes,
                          args, expr_available_ids, gtf_tcons_to_loc,
                          gwas_pvalues, term_matches, query_info, go_map,
                          all_rising):
    """Run the shared downstream pipeline for one expanded module.

    Loads expression, builds the co-expression network, computes centrality,
    identifies rising hubs, and writes output files.  Appends rising-hub
    results to ``all_rising``.
    """
    if len(expanded) < 2:
        print(f"  [{mid}] SKIP: insufficient expanded isoforms (<2)")
        return

    # Load expression data for expanded set (isoform-level)
    expr_df_iso, expr_id_map = load_expression_subset(
        args.expression,
        expanded,
    )

    if expr_df_iso.empty or len(expr_df_iso) < 2:
        print(f"  [{mid}] SKIP: insufficient expression data (<2 isoforms)")
        return

    # Gene-completeness: for genes that have at least one isoform in
    # the expanded set, also load any missing siblings.
    existing_keys = set(expr_df_iso.index)
    locs_in_expr = {k.split(':')[0] for k in existing_keys if ':' in k}
    missing_isoforms = set()
    for loc_id in locs_in_expr:
        for ek in expr_available_ids:
            if ':' in ek and ek.startswith(loc_id + ':') and ek not in existing_keys:
                missing_isoforms.add(ek.split(':')[1])
    if missing_isoforms:
        expr_df_sibs, sib_id_map = load_expression_subset(
            args.expression, missing_isoforms)
        if not expr_df_sibs.empty:
            new_keys = set(expr_df_sibs.index) - existing_keys
            expr_df_sibs = expr_df_sibs.loc[list(new_keys)]
            if not expr_df_sibs.empty:
                expr_df_iso = pd.concat([expr_df_iso, expr_df_sibs])
                for qid, keys in sib_id_map.items():
                    expr_id_map.setdefault(qid, []).extend(keys)
                print(f"  [EXPR]   added {len(expr_df_sibs)} missing sibling "
                      f"isoforms for {len(locs_in_expr)} genes")

    # Aggregate to gene-level expression
    expr_df_gene, gene_to_isoforms = aggregate_to_gene_level(expr_df_iso)

    if expr_df_gene.empty or len(expr_df_gene) < 2:
        print(f"  [{mid}] SKIP: insufficient genes after aggregation (<2)")
        return

    # Build co-expression network on gene-level expression
    G, expr_filt = build_coexpression_network(
        expr_df_gene,
        corr_threshold=args.corr_threshold,
        fdr_threshold=args.fdr_threshold,
        min_mean_tpm=args.min_tpm,
    )

    if G.number_of_nodes() < 2:
        print(f"  [{mid}] SKIP: insufficient nodes in co-expression network")
        return

    # Force-connect: for force-included genes that ended up isolated
    if args.force_genes:
        forced_set = {g.strip() for g in args.force_genes.split(',') if g.strip()}
        for fg in forced_set & set(G.nodes()):
            if G.degree(fg) > 0:
                continue
            fg_idx = expr_filt.index.get_loc(fg)
            fg_vec = expr_filt.iloc[fg_idx].values.astype(np.float64)
            best_r, best_partner = 0.0, None
            for other in G.nodes():
                if other == fg: continue
                o_vec = expr_filt.iloc[expr_filt.index.get_loc(other)].values.astype(np.float64)
                r = np.corrcoef(fg_vec, o_vec)[0, 1]
                if abs(r) > abs(best_r):
                    best_r, best_partner = r, other
            if best_partner and abs(best_r) >= 0.4:
                G.add_edge(fg, best_partner, weight=best_r, rho=best_r)
                print(f"  [FORCE] {fg} connected to {best_partner} "
                      f"(r={best_r:.3f}, relaxed from "
                      f"global |ρ|≥{args.corr_threshold})")
            else:
                print(f"  [FORCE] {fg} could not be connected "
                      f"(best |r|={abs(best_r):.3f} < 0.4)")

    # Multi-hop expansion
    if args.expansion_hops >= 2 and args.expansion_link_rho > 0:
        orig_module_genes = set()
        for iso in mod_isos_in_nodes:
            for ek in expr_id_map.get(iso, []):
                if ':' in ek:
                    orig_module_genes.add(ek.split(':')[0])
                else:
                    orig_module_genes.add(ek)
            if iso.startswith('LOC'):
                orig_module_genes.add(iso)

        seed_genes = orig_module_genes & set(G.nodes())
        n_seeds = len(seed_genes)
        print(f"  [HOPS] hop 1 → {n_seeds} original module genes in network")

        total_added = 0
        current_seeds = set(seed_genes)
        all_linked = set()

        for hop in range(2, args.expansion_hops + 1):
            linked = expand_by_coexpression_links(
                G, current_seeds,
                link_rho=args.expansion_link_rho,
                max_neighbors=args.max_expansion,
            )
            new_to_add = linked - set(expr_df_gene.index) - all_linked
            if not new_to_add:
                print(f"  [HOPS] hop {hop} → no new genes found "
                      f"(all {len(linked)} linked genes already present)")
                break

            all_linked |= new_to_add
            current_seeds = new_to_add
            total_added += len(new_to_add)
            print(f"  [HOPS] hop {hop} → {len(new_to_add)} new genes "
                  f"via co-expression (|rho| >= {args.expansion_link_rho})")

        if all_linked:
            linked_query_ids = set()
            for gid in all_linked:
                linked_query_ids.add(gid)
                for ek in gene_to_isoforms.get(gid, []):
                    if ':' in ek:
                        linked_query_ids.add(ek.split(':')[1])

            expr_df_linked, linked_id_map = load_expression_subset(
                args.expression,
                linked_query_ids,
            )

            if not expr_df_linked.empty:
                expr_df_iso = pd.concat([expr_df_iso, expr_df_linked])
                expr_df_gene, gene_to_isoforms = aggregate_to_gene_level(
                    expr_df_iso)
                print(f"  [COEXPR] rebuilding network with "
                      f"{len(expr_df_gene)} genes ({total_added} added via hops) →")
                G, expr_filt = build_coexpression_network(
                    expr_df_gene,
                    corr_threshold=args.corr_threshold,
                    fdr_threshold=args.fdr_threshold,
                    min_mean_tpm=args.min_tpm,
                )
                if args.force_genes:
                    forced_set = {g.strip() for g in args.force_genes.split(',') if g.strip()}
                    for fg in forced_set & set(G.nodes()):
                        if G.degree(fg) > 0: continue
                        fg_vec = expr_filt.iloc[expr_filt.index.get_loc(fg)].values.astype(np.float64)
                        best_r, best_partner = 0.0, None
                        for other in G.nodes():
                            if other == fg: continue
                            o_vec = expr_filt.iloc[expr_filt.index.get_loc(other)].values.astype(np.float64)
                            r = np.corrcoef(fg_vec, o_vec)[0, 1]
                            if abs(r) > abs(best_r):
                                best_r, best_partner = r, other
                        if best_partner and abs(best_r) >= 0.4:
                            G.add_edge(fg, best_partner, weight=best_r, rho=best_r)
                            print(f"  [FORCE] (hops) {fg} connected to {best_partner} "
                                  f"(r={best_r:.3f})")

                for qid, keys in linked_id_map.items():
                    expr_id_map.setdefault(qid, []).extend(keys)
            else:
                print(f"  [HOPS] no expression data for linked genes")
    else:
        print(f"  [HOPS] single-hop (annotation only, "
              f"use --expansion-hops ≥2 to follow co-expression links)")

    if G.number_of_nodes() < 2:
        print(f"  [{mid}] SKIP: insufficient nodes after multi-hop expansion")
        return

    # Compute centrality metrics (gene-level)
    print(f"  [CENTR] Computing centrality metrics →")
    centrality_df = compute_centrality_metrics(G, expr_filt)

    # Build gene-ID → best GWAS p-value mapping
    gene_to_gwas_p = {}
    for gene_id in centrality_df.index:
        best_p = 1.0
        isoforms = gene_to_isoforms.get(gene_id, [gene_id])
        for iso_key in isoforms:
            for qid, mapped_keys in expr_id_map.items():
                if iso_key in mapped_keys:
                    p = gwas_pvalues.get(qid, 1.0)
                    if p < best_p:
                        best_p = p
                    if qid.startswith('LOC') and qid in gwas_pvalues:
                        p = gwas_pvalues[qid]
                        if p < best_p:
                            best_p = p
        p = gwas_pvalues.get(gene_id, best_p)
        if p < best_p:
            best_p = p
        gene_to_gwas_p[gene_id] = best_p

    # Build gene-level expression-key → query-ID mapping
    gene_key_to_query = {}
    for qid, expr_keys in expr_id_map.items():
        for ek in expr_keys:
            if ':' in ek:
                gene = ek.split(':', 1)[0]
            else:
                gene = ek
            existing = gene_key_to_query.get(gene, '')
            if not existing or (qid.startswith('TCONS') and not existing.startswith('TCONS')):
                gene_key_to_query[gene] = qid

    # Compute combined regulatory scores (gene-level)
    print(f"  [SCORE] Computing combined regulatory scores →")
    ranking_df = compute_regulatory_scores(
        centrality_df, gene_to_gwas_p, gene_key_to_query,
        gwas_weight=args.gwas_weight,
        expr_weight=args.expr_weight,
        expression_only=args.expression_only,
    )

    # Identify rising hubs
    rising_df = identify_rising_hubs(
        ranking_df,
        min_rank_rise=args.min_rank_rise,
        top_n=30,
        min_gwas_pvalue=args.min_gwas_pvalue,
    )
    print(f"  [RISE] {len(rising_df)} rising hubs detected")
    if not rising_df.empty:
        print(f"    Top 5:")
        for iso, row in rising_df.head(5).iterrows():
            print(f"      {iso:40s}  Δrank=+{int(row['rank_change']):3d}  "
                  f"p={row['gwas_pvalue']:.2e}  kME={row['kME']:.3f}")

    # Write outputs
    write_module_outputs(
        mid, mod_name, expanded, ranking_df, rising_df,
        G, expr_filt, term_matches, query_info,
        gene_to_isoforms, gene_key_to_query, args.out_prefix,
    )

    if not rising_df.empty:
        all_rising.append((mid, mod_name, rising_df))


def extract_module_terms(module_isoforms, nodes, modules_df_row=None):
    """Extract characterizing functional terms for a module.

    Aggregates PFAM, GO, KEGG, and COG terms from module members and ranks
    them by prevalence (how many isoforms in the module have each term).
    Filters out generic terms.

    Returns
    -------
    dict
        Keys: 'pfam', 'go', 'kegg', 'cog', 'desc' ? list of (term, count) tuples
        sorted by count descending.  Only includes terms shared by ?2 isoforms.
    """
    pfam_counter = Counter()
    go_counter = Counter()
    kegg_counter = Counter()
    cog_counter = Counter()
    desc_counter = Counter()

    for iso in module_isoforms:
        nd = nodes.get(iso, {})

        # Helper: safely get a string value (handle NaN)
        def _str(val):
            if val is None or (isinstance(val, float) and np.isnan(val)):
                return ''
            return str(val)

        # PFAM
        pfam_text = _str(nd.get('PFAM', ''))
        if pfam_text:
            for t in pfam_text.split(','):
                t = t.strip()
                if t and t != '-':
                    pfam_counter[t] += 1

        # GO descriptions (already resolved in nodes)
        go_text = _str(nd.get('GO_description', ''))
        if go_text:
            for t in go_text.split(';'):
                t = t.strip()
                if t and t != '-' and not _is_generic_go(t):
                    go_counter[t] += 1

        # KEGG
        kegg_text = _str(nd.get('KEGG_Pathway', ''))
        if kegg_text:
            for t in kegg_text.split(';'):
                t = t.strip()
                if t and t != '-' and not _is_broad_kegg(t):
                    kegg_counter[t] += 1

        # COG
        cog_text = _str(nd.get('COG_category', ''))
        if cog_text:
            for t in cog_text.split(';'):
                t = t.strip()
                if t and t != '-' and t not in _GENERIC_COG_EDGE_TERMS:
                    cog_counter[t] += 1

        # Description
        desc_text = _str(nd.get('Annotation', ''))
        if desc_text and len(desc_text) > 3 and not _is_vague_description(desc_text):
            desc_key = desc_text[:120].strip().lower()
            if desc_key:
                desc_counter[desc_key] += 1

    # Keep terms shared by ≥2 isoforms in the module.
    # For TF modules we relax PFAM & DESC to ≥1 so highly specific terms
    # (Myb_CC_LHEQLE, GLK1-like) that appear in only one seed are kept.
    # Never relax DESC for non-TF modules — substring matching on singleton
    # phrases like "repeat-containing protein" inflates scores across
    # thousands of unrelated genes and crowds out genuine module members.
    n_seeds = len(module_isoforms)
    pfam_terms_module = {t for t in pfam_counter}
    is_tf_mod = bool(pfam_terms_module & _TF_PFAM_TERMS)
    pfam_min = 1 if (n_seeds <= 5 and is_tf_mod) else 2
    desc_min = 1 if (n_seeds <= 5 and is_tf_mod) else 2
    result = {
        'pfam': [(t, c) for t, c in pfam_counter.most_common(20) if c >= pfam_min],
        'go': [(t, c) for t, c in go_counter.most_common(20) if c >= 2],
        'kegg': [(t, c) for t, c in kegg_counter.most_common(10) if c >= 2],
        'cog': [(t, c) for t, c in cog_counter.most_common(10) if c >= 2],
        'desc': [(t, c) for t, c in desc_counter.most_common(10) if c >= desc_min],
    }
    return result


# ---------------------------------------------------------------------------
# Module expansion via annotation search
# ---------------------------------------------------------------------------

def _scan_expression_ids(expression_path):
    """Quickly scan the expression matrix for all available gene/isoform IDs.

    Returns
    -------
    expr_ids : set of str
        All LOC:TCONS keys, plus extracted TCONS and LOC parts.
    """
    expr_ids = set()
    with open(expression_path, 'r') as fh:
        _ = fh.readline()  # skip header
        for line in fh:
            key = line.split('\t', 1)[0]
            expr_ids.add(key)
            if ':' in key:
                loc, tcons = key.split(':', 1)
                expr_ids.add(loc)
                expr_ids.add(tcons)
    return expr_ids


def expand_module(module_terms, annotation_index, query_info,
                  go_map=None, min_shared_terms=2, max_expansion=500,
                  expr_available_ids=None):
    """Expand a module by searching the full annotation for functionally
    related isoforms.

    For each characterizing term, retrieves the set of query_ids from the
    annotation index.  An isoform must share ? ``min_shared_terms`` with the
    module's term set to be included.

    PFAM terms are weighted most heavily (they are the most specific),
    followed by specific GO terms, then KEGG pathway membership.

    Isoforms present in the expression dataset receive a bonus so they
    survive the max_expansion cutoff : this prevents expression-relevant
    genes from being crowded out by annotation-only entries.

    Parameters
    ----------
    module_terms : dict
        Output of ``extract_module_terms``.
    annotation_index : dict
        Output of ``build_annotation_index``.
    query_info : dict
        query_id ? annotation info dict.
    go_map : dict or None
        GO ID ? description mapping (to resolve GO IDs to text for matching).
    min_shared_terms : int
        Minimum shared terms for inclusion.
    max_expansion : int
        Maximum expanded module size (top-N by shared term count).
    expr_available_ids : set or None
        Pre-scanned set of gene/isoform IDs present in the expression matrix.
        Used as a tie-breaking bonus so expression-relevant genes are kept.

    Returns
    -------
    expanded : set of str
        Query IDs (TCONS_* or LOC*) in the expanded module.
    term_matches : dict
        query_id ? list of matched term strings.
    """
    # Collect the actual term strings to match against
    pfam_terms = {t for t, _ in module_terms.get('pfam', [])}
    go_desc_terms = {t for t, _ in module_terms.get('go', [])}
    kegg_terms = {t for t, _ in module_terms.get('kegg', [])}
    desc_terms = {t for t, _ in module_terms.get('desc', [])}

    # Collect candidate query_ids from each term's index entries
    candidates = set()
    for term in pfam_terms:
        candidates |= annotation_index['pfam'].get(term, set())
    for term in desc_terms:
        candidates |= annotation_index['desc'].get(term, set())

    # For GO terms, we need to match by GO ID in the annotation
    # The module_terms['go'] contains human-readable descriptions.
    # We need to reverse-lookup GO IDs from descriptions using go_map.
    if go_map:
        desc_to_id = {v.lower(): k for k, v in go_map.items()}
        for go_desc in go_desc_terms:
            go_id = desc_to_id.get(go_desc.lower())
            if go_id:
                candidates |= annotation_index['go'].get(go_id, set())

    # Score ALL candidates.  Two-tier sort (score then expr-flag) ensures
    # expression-available genes ALWAYS rank above annotation-only genes
    # at the same score.  For modules with very common terms (e.g. PPR
    # with 3,000+ expression-available candidates), we also pre-filter
    # to expression-available only : the whole point is expression-based
    # analysis.  A fallback prevents pathological cases where too few
    # expression-available candidates exist.
    if expr_available_ids:
        candidates_expr = candidates & expr_available_ids
        if len(candidates_expr) >= 50:
            candidates = candidates_expr
        else:
            print(f"  [EXPAND] only {len(candidates_expr)} expression-available "
                  f"candidates : falling back to full annotation pool "
                  f"({len(candidates)} total)")

    scored = []
    for qid in candidates:
        info = query_info.get(qid, {})
        score = 0.0
        matched = []

        # PFAM match
        for t in info.get('pfam', []):
            if t in pfam_terms:
                score += 2  # PFAM is strong signal
                matched.append(f"PFAM:{t}")

        # GO ID match (need to resolve via go_map)
        for go_id in info.get('go', []):
            if go_map and go_id in go_map:
                go_desc = go_map[go_id].lower()
                for mod_go in go_desc_terms:
                    if mod_go.lower() in go_desc or go_desc in mod_go.lower():
                        score += 1
                        matched.append(f"GO:{go_desc}")
                        break

        # Description fuzzy match
        desc = info.get('desc', '')
        if desc:
            for mod_desc in desc_terms:
                if mod_desc[:60] in desc.lower() or desc.lower()[:60] in mod_desc:
                    score += 1
                    matched.append(f"DESC:{mod_desc[:60]}")
                    break

        if score >= min_shared_terms:
            scored.append((qid, score, matched))

    # Sort: primary = annotation score DESC,
    #        secondary = matched term count DESC (more specific matches),
    #        tertiary = expression-available (True > False),
    #        quaternary = query_id (deterministic tiebreaker).
    # This prevents random exclusion when thousands of candidates tie
    # at the same score (e.g. all 3,134 PPR genes score 7).
    scored.sort(key=lambda x: (-x[1], -len(x[2]), not (expr_available_ids and x[0] in expr_available_ids), x[0]))

    top_n = scored[:max_expansion]

    # Rescue: expression-available isoforms that scored well but fell outside
    # the top max_expansion (because two-tier sort already prioritized
    # expr-available genes, this catches highly-scoring genes from narrow
    # annotation terms that share fewer module-characterizing terms but are
    # still functionally related).
    rescue_cap = max_expansion
    rescued = []
    rescue_threshold = max(min_shared_terms + 1, 3)
    for qid, score, matched in scored[max_expansion:]:
        if len(rescued) >= rescue_cap:
            break
        if expr_available_ids and qid in expr_available_ids:
            if score >= rescue_threshold:
                rescued.append((qid, score, matched))

    if rescued:
        print(f"  [EXPAND] rescued {len(rescued)} expression-available isoforms "
              f"that scored ?{rescue_threshold} but fell outside top {max_expansion} "
              f"(capped at {rescue_cap})")

    scored = top_n + rescued
    expanded = {qid for qid, _, _ in scored}
    term_matches = {qid: matches for qid, _, matches in scored}

    return expanded, term_matches


# ---------------------------------------------------------------------------
# Expression data loading (isoform-level and gene-level)
# ---------------------------------------------------------------------------

def load_expression_subset(expression_path, target_ids, gtf_tcons_to_loc=None):
    """Load expression data for a target set of isoforms (isoform-level).

    Returns every matching transcript row : use `aggregate_to_gene_level`
    afterwards to collapse isoforms into gene-level TPMs for co-expression.

    The expression file has format:
        geneID:transcriptID<TAB>TPM1<TAB>TPM2<...>

    Returns
    -------
    expr_df : pd.DataFrame
        Rows = isoform expression keys (LOC:TCONS), columns = sample names.
    id_mapping : dict
        Maps annotation query_id ? list of expression row keys.
    """
    print(f"  [EXPR] Loading expression data for {len(target_ids)} target isoforms ?")
    t0 = time.time()

    tcons_targets = {t for t in target_ids if t.startswith('TCONS')}
    loc_targets = {t for t in target_ids if t.startswith('LOC')}

    rows = []
    id_mapping = defaultdict(list)  # query_id ? list of expression row keys

    with open(expression_path, 'r') as fh:
        header = fh.readline().rstrip('\n').split('\t')
        sample_names = header[1:]

        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 2:
                continue
            expr_key = parts[0]  # e.g., "LOC111886547:TCONS_00184145"

            if ':' in expr_key:
                loc_part, tcons_part = expr_key.split(':', 1)
            else:
                loc_part, tcons_part = expr_key, expr_key

            matched = False
            if tcons_part in tcons_targets:
                matched = True
                id_mapping[tcons_part].append(expr_key)
            if loc_part in loc_targets:
                matched = True
                id_mapping[loc_part].append(expr_key)

            if matched:
                try:
                    values = [float(v) for v in parts[1:]]
                except ValueError:
                    values = [0.0] * len(sample_names)
                rows.append([expr_key] + values)

    expr_df = pd.DataFrame(rows, columns=['geneID:transcriptID'] + list(sample_names))
    if not expr_df.empty:
        expr_df = expr_df.set_index('geneID:transcriptID')

    elapsed = time.time() - t0
    n_found = len(expr_df)
    n_mapped = len(id_mapping)
    print(f"  [EXPR]   {n_found} expression rows for {n_mapped} query IDs "
          f"({elapsed:.1f}s)")

    return expr_df, dict(id_mapping)


def aggregate_to_gene_level(expr_df):
    """Sum all transcript-level TPMs into gene-level TPMs.

    Skips XLOC_* genes (StringTie-only intergenic predictions, no RefSeq
    mapping). Input: ``LOC111886547:TCONS_00184145``. Output: ``LOC111886547``
    (one row per gene, TPMs summed across all transcripts).

    Returns
    -------
    gene_expr : pd.DataFrame
        Rows = gene IDs (LOC*), columns = sample names.
    gene_to_transcripts : dict
        gene_id --> list of isoform expression keys that contributed to it.
    """
    if expr_df.empty:
        return expr_df, {}

    gene_to_transcripts = defaultdict(list)
    gene_sums = defaultdict(lambda: None)

    n_skipped = 0
    for idx in expr_df.index:
        if ':' in idx:
            loc = idx.split(':', 1)[0]
        else:
            loc = idx
        # Skip XLOC_*: StringTie-only intergenic clusters, not real genes
        if loc.startswith('XLOC_'):
            n_skipped += 1
            continue
        gene_to_transcripts[loc].append(idx)

        vals = expr_df.loc[idx].values.astype(np.float64)
        if gene_sums[loc] is None:
            gene_sums[loc] = vals.copy()
        else:
            gene_sums[loc] += vals

    gene_expr = pd.DataFrame(
        {loc: gene_sums[loc] for loc in gene_sums},
        index=expr_df.columns,
    ).T
    gene_expr.index.name = 'gene_id'

    n_isoforms = len(expr_df)
    n_genes = len(gene_expr)
    suff = f"  ({n_skipped} XLOC skipped)" if n_skipped > 0 else ""
    print(f"  [EXPR]   aggregated {n_isoforms} isoforms ? {n_genes} genes "
          f"({n_isoforms - n_genes} isoforms collapsed){suff}")

    return gene_expr, dict(gene_to_transcripts)


# ---------------------------------------------------------------------------
# Co-expression network construction
# ---------------------------------------------------------------------------

def build_coexpression_network(expr_df, corr_threshold=0.7, fdr_threshold=0.05,
                                min_mean_tpm=0.1, max_genes=1500):
    """Build a co-expression network from a TPM dataframe.

    1. Filters genes with mean TPM < min_mean_tpm
    2. If > max_genes remain, subsamples to top max_genes by mean TPM
       (highest-expressed genes are most informative for co-expression)
    3. Computes all pairwise Spearman rank correlations
    4. Applies Benjamini-Hochberg FDR correction
    5. Returns edges passing both |?| > corr_threshold and FDR < fdr_threshold

    Parameters
    ----------
    expr_df : pd.DataFrame
        Rows = genes, columns = samples (TPM values).
    corr_threshold : float
        Minimum absolute Spearman ? for edge inclusion.
    fdr_threshold : float
        Maximum FDR-adjusted p-value for edge inclusion.
    min_mean_tpm : float
        Minimum mean TPM across samples to retain a gene.
    max_genes : int
        Maximum number of genes to carry forward to pairwise correlation.
        Beyond this, O(n?) becomes intractable.  Only the highest-expressed
        genes are kept : low-expression noise is filtered out.

    Returns
    -------
    G : nx.Graph
        Co-expression graph (nodes = gene IDs, edges with weight=?).
    expr_filtered : pd.DataFrame
        Filtered expression data (only retained genes).
    """
    # Filter low-expression isoforms
    mean_tpm = expr_df.mean(axis=1)
    keep = mean_tpm >= min_mean_tpm
    expr_filt = expr_df.loc[keep]

    # Hard cap: if too many genes, subsample to top by mean TPM
    # so pairwise correlations stay tractable (O(n?) cost)
    if len(expr_filt) > max_genes:
        top_idx = mean_tpm[keep].sort_values(ascending=False).head(max_genes).index
        expr_filt = expr_filt.loc[top_idx]
        print(f"  [COEXPR]   subsampled to top {max_genes} genes by mean TPM "
              f"(from {keep.sum()}) for tractable O(n?) correlation")

    if len(expr_filt) < 2:
        print(f"  [COEXPR]   WARNING: only {len(expr_filt)} isoforms after "
              f"TPM ? {min_mean_tpm} filter : cannot build network")
        return nx.Graph(), expr_filt

    n_genes = len(expr_filt)
    n_pairs = n_genes * (n_genes - 1) // 2
    print(f"  [COEXPR]   {n_genes} isoforms after TPM ? {min_mean_tpm} filter "
          f"? {n_pairs:,} pairwise correlations")

    # Compute pairwise Spearman correlations
    gene_ids = expr_filt.index.tolist()
    expr_matrix = expr_filt.values  # (n_genes, n_samples)

    edges = []
    pvalues = []
    edge_records = []

    # Use scipy's spearmanr in a loop
    for i in range(n_genes):
        for j in range(i + 1, n_genes):
            rho, pval = spearmanr(expr_matrix[i], expr_matrix[j])
            if not np.isnan(rho):
                edge_records.append({
                    'source': gene_ids[i],
                    'target': gene_ids[j],
                    'rho': rho,
                    'pvalue': pval,
                })
                pvalues.append(pval)

    if not edge_records:
        print(f"  [COEXPR]   WARNING: no valid correlations computed")
        return nx.Graph(), expr_filt

    # Benjamini-Hochberg FDR correction
    pvalues = np.array(pvalues)
    n_tests = len(pvalues)
    sorted_indices = np.argsort(pvalues)
    sorted_p = pvalues[sorted_indices]
    ranks = np.arange(1, n_tests + 1)
    bh_critical = ranks / n_tests * fdr_threshold

    # Find the largest p-value that is <= its BH critical value
    significant = sorted_p <= bh_critical
    if significant.any():
        max_sig_idx = np.where(significant)[0][-1]
        pvalue_cutoff = sorted_p[max_sig_idx]
    else:
        pvalue_cutoff = 0.0

    print(f"  [COEXPR]   BH FDR threshold: p ? {pvalue_cutoff:.2e} "
          f"(FDR = {fdr_threshold})")

    # Build graph with edges passing both thresholds
    G = nx.Graph()
    for iso_id in gene_ids:
        G.add_node(iso_id)

    n_sig = 0
    for rec in edge_records:
        if abs(rec['rho']) >= corr_threshold and rec['pvalue'] <= pvalue_cutoff:
            G.add_edge(rec['source'], rec['target'],
                       weight=rec['rho'], pvalue=rec['pvalue'])
            n_sig += 1

    print(f"  [COEXPR]   {n_sig} significant edges "
          f"(|rho| >= {corr_threshold}, FDR < {fdr_threshold})")

    return G, expr_filt


def expand_by_coexpression_links(G, seed_genes, link_rho=0.6, max_neighbors=500):
    """Expand a gene set by following strong co-expression edges from seed genes.

    For each hop, adds all neighbors connected to any seed gene with
    |rho| >= link_rho.  This captures genes that were not in the original
    annotation-based expansion but are tightly co-expressed with module members.

    Parameters
    ----------
    G : nx.Graph
        Co-expression graph (nodes = gene IDs, edge weight = Spearman rho).
    seed_genes : set of str
        Genes to expand from (original module members + annotation expansion).
    link_rho : float
        Minimum |rho| for a co-expression edge to trigger inclusion.
    max_neighbors : int
        Maximum number of neighbor genes to add (prevents runaway expansion).

    Returns
    -------
    new_genes : set of str
        Genes that should be added to the module (not already in seed_genes).
    """
    new_genes = set()
    for seed in seed_genes:
        if seed not in G:
            continue
        for neighbor in G.neighbors(seed):
            if neighbor in seed_genes:
                continue
            edge_data = G.get_edge_data(seed, neighbor)
            rho = abs(edge_data.get('weight', 0))
            if rho >= link_rho:
                new_genes.add(neighbor)
        if len(new_genes) >= max_neighbors:
            break

    return new_genes


# ---------------------------------------------------------------------------
# Centrality metrics
# ---------------------------------------------------------------------------

def compute_centrality_metrics(G, expr_df):
    """Compute multiple centrality metrics for a co-expression graph.

    Parameters
    ----------
    G : nx.Graph
        Co-expression graph.
    expr_df : pd.DataFrame
        Filtered expression data (rows = isoforms, cols = samples).

    Returns
    -------
    metrics : pd.DataFrame
        Index = isoform ID, columns:
        - degree: number of co-expression edges
        - betweenness: betweenness centrality (normalized)
        - eigenvector: eigenvector centrality
        - clustering: local clustering coefficient
        - kME: module eigengene connectivity (correlation to PC1)
    """
    if G.number_of_nodes() == 0:
        return pd.DataFrame()

    nodes = list(G.nodes())

    # Degree centrality (normalized by n-1)
    degree = {n: G.degree(n) / max(G.number_of_nodes() - 1, 1) for n in nodes}

    # Betweenness centrality (normalized)
    betweenness = nx.betweenness_centrality(G, normalized=True)

    # Eigenvector centrality (may fail on disconnected graphs)
    try:
        eigenvector = nx.eigenvector_centrality_numpy(G, max_iter=500)
    except (nx.PowerIterationFailedConvergence, RuntimeError, TypeError):
        # Fallback 1 (small/disconnected graphs): power-iteration method
        try:
            eigenvector = nx.eigenvector_centrality(G, max_iter=500)
        except (nx.PowerIterationFailedConvergence, RuntimeError):
            # Fallback 2: use weighted degree
            eigenvector = {n: G.degree(n, weight='weight') / max(
                sum(d.get('weight', 1) for _, _, d in G.edges(data=True)), 1)
                for n in nodes}

    # Clustering coefficient
    clustering = nx.clustering(G, weight='weight')

    # kME: module eigengene connectivity
    # Compute first principal component of expression data via SVD (numpy,
    # no sklearn dependency).
    kME = {}
    if len(expr_df) >= 2 and expr_df.shape[1] >= 2:
        expr_subset = expr_df.loc[expr_df.index.isin(nodes)]
        if len(expr_subset) >= 2:
            try:
                # expr_subset: rows=genes, cols=samples
                # Center the data
                X = expr_subset.values.astype(np.float64).T  # (samples, genes)
                X -= X.mean(axis=0)
                # SVD for PC1
                U, S, Vt = np.linalg.svd(X, full_matrices=False)
                pc1 = U[:, 0]  # first left singular vector = sample scores on PC1
                for i, iso_id in enumerate(expr_subset.index):
                    gene_vec = expr_subset.iloc[i].values.astype(np.float64)
                    if np.std(gene_vec) > 0:
                        r, _ = spearmanr(gene_vec, pc1)
                        kME[iso_id] = abs(r) if not np.isnan(r) else 0.0
                    else:
                        kME[iso_id] = 0.0
            except np.linalg.LinAlgError:
                for iso_id in nodes:
                    kME[iso_id] = 0.0
        else:
            for iso_id in nodes:
                kME[iso_id] = 0.0
    else:
        for iso_id in nodes:
            kME[iso_id] = 0.0

    # Build dataframe
    metrics = pd.DataFrame({
        'degree_centrality': degree,
        'betweenness_centrality': betweenness,
        'eigenvector_centrality': eigenvector,
        'clustering_coefficient': clustering,
        'kME': kME,
    })
    metrics.index.name = 'isoform'
    return metrics


# ---------------------------------------------------------------------------
# GWAS p-value lookup
# ---------------------------------------------------------------------------

def load_gwas_pvalues(qqman_path, gtf_tcons_to_loc=None):
    """Load GWAS p-values from a QQMan file.

    Returns
    -------
    dict
        Maps gene_id (TCONS or LOC) ? best p-value.
    """
    print(f"[GWAS] Loading p-values from: {qqman_path}")
    df = pd.read_csv(
        qqman_path, sep='\t', header=None,
        names=['original_id', 'chr', 'pos', 'pvalue'],
        dtype={'original_id': str, 'chr': str, 'pos': int, 'pvalue': float}
    )

    # Extract canonical gene ID
    re_loc = re.compile(r'LOC\d+')
    re_tcons = re.compile(r'TCONS_\d+')

    def _extract(raw):
        m = re_loc.search(raw)
        if m:
            return m.group(0)
        m = re_tcons.search(raw)
        if m:
            return m.group(0)
        return raw

    df['gene_id'] = df['original_id'].apply(_extract)

    # Keep best p-value per gene_id
    best_p = df.groupby('gene_id')['pvalue'].min().to_dict()
    print(f"[GWAS]   {len(best_p):,} unique gene IDs with p-values")
    return best_p


# ---------------------------------------------------------------------------
# Combined regulatory score
# ---------------------------------------------------------------------------

def compute_regulatory_scores(centrality_df, gwas_pvalues, expr_key_to_query,
                               gwas_weight=0.3, expr_weight=0.7,
                               expression_only=False):
    """Compute combined regulatory score integrating GWAS + expression evidence.

    For each isoform:
        Z_score = rho (kME)

    Missing GWAS p-values are assigned p=1.0 (no evidence).

    Parameters
    ----------
    centrality_df : pd.DataFrame
        Output of ``compute_centrality_metrics``.
    gwas_pvalues : dict
        gene_id ? p-value mapping.
    expr_key_to_query : dict
        Maps expression key (LOC:TCONS) ? annotation query_id.
    gwas_weight : float
        Weight for GWAS signal in combined score.
    expr_weight : float
        Weight for expression metrics (split between betweenness and kME).

    Returns
    -------
    ranking_df : pd.DataFrame
        Sorted by combined_score descending.
    """
    if centrality_df.empty:
        return centrality_df

    df = centrality_df.copy()

    # Map expression keys to query IDs and look up GWAS p-values
    # At gene level: expression keys are LOC IDs, gwas_pvalues are already
    # keyed by gene ID from gene_to_gwas_p.
    def _get_pval(expr_key):
        # Direct lookup in gwas_pvalues (gene-level)
        pval = gwas_pvalues.get(expr_key)
        if pval is not None:
            return pval
        # Try via expr_key_to_query mapping
        qid = expr_key_to_query.get(expr_key, expr_key)
        pval = gwas_pvalues.get(qid)
        if pval is not None:
            return pval
        return 1.0

    df['gwas_pvalue'] = [_get_pval(k) for k in df.index]
    # Guard against p=0
    df['gwas_pvalue'] = df['gwas_pvalue'].clip(lower=1e-300)
    df['neg_log10_p'] = -np.log10(df['gwas_pvalue'])

    # Evidence type: does this gene have any GWAS signal?
    # p=1.0 or p=1e-300 means no GWAS hit was found for this gene.
    df['has_gwas_signal'] = df['gwas_pvalue'] < 1.0

    # Z-score normalize each metric
    def _zscore(s):
        std = s.std()
        if std == 0:
            return pd.Series(0.0, index=s.index)
        return (s - s.mean()) / std

    z_neg_log_p = _zscore(df['neg_log10_p'])
    z_betweenness = _zscore(df['betweenness_centrality'])
    z_kme = _zscore(df['kME'])

    # In expression-only mode, score purely by co-expression centrality.
    # GWAS contributes nothing : this surfaces genes that are GWAS-invisible
    # (essential genes, redundant regulators, post-transcriptional factors).
    if expression_only:
        alpha = 0.0
        beta = 0.6   # betweenness
        gamma = 0.4  # kME
    else:
        alpha = gwas_weight
        beta = expr_weight * 0.6
        gamma = expr_weight * 0.4

    df['z_gwas'] = z_neg_log_p
    df['z_betweenness'] = z_betweenness
    df['z_kME'] = z_kme
    df['combined_score'] = (alpha * z_neg_log_p +
                            beta * z_betweenness +
                            gamma * z_kme)

    # Evidence label
    df['evidence'] = df['has_gwas_signal'].apply(
        lambda x: 'gwas+expr' if x else 'expr_only')

    # Rank by combined score
    df = df.sort_values('combined_score', ascending=False)
    df['rank_combined'] = range(1, len(df) + 1)

    # Rank by GWAS alone for comparison
    df = df.sort_values('neg_log10_p', ascending=False)
    df['rank_gwas'] = range(1, len(df) + 1)

    # Sort back by combined score
    df = df.sort_values('combined_score', ascending=False)

    # Rank change: positive = rose in combined rank vs GWAS rank
    # In expression-only mode (GWAS weight=0), this measures how much
    # expression centrality lifts a gene above its GWAS baseline.
    df['rank_change'] = df['rank_gwas'] - df['rank_combined']

    return df


# ---------------------------------------------------------------------------
# Rising hub detection
# ---------------------------------------------------------------------------

def identify_rising_hubs(ranking_df, min_rank_rise=10, top_n=30,
                         min_gwas_pvalue=0.01):
    """Identify "rising hubs" : isoforms that rank much higher by combined
    score than by GWAS alone.

    These are genes that GWAS would have missed or deprioritized, but
    expression network analysis flags as central regulators.

    Parameters
    ----------
    ranking_df : pd.DataFrame
        Output of ``compute_regulatory_scores``.
    min_rank_rise : int
        Minimum rank improvement to be considered a rising hub.
    top_n : int
        Maximum number of rising hubs to return.
    min_gwas_pvalue : float
        Maximum GWAS p-value to be a "rising hub" (must be GWAS-weak).

    Returns
    -------
    pd.DataFrame
        Subset of ranking_df sorted by rank_change descending.
    """
    rising = ranking_df[
        (ranking_df['rank_change'] >= min_rank_rise) &
        (ranking_df['gwas_pvalue'] >= min_gwas_pvalue)
    ].copy()
    rising = rising.sort_values('rank_change', ascending=False).head(top_n)
    return rising


# ---------------------------------------------------------------------------
# Module deduplication (consolidate multi-isoform genes)
# ---------------------------------------------------------------------------

def _extract_tcons_from_expr_key(expr_key):
    """Extract TCONS part from expression key 'LOC:TCONS'."""
    if ':' in expr_key:
        return expr_key.split(':')[1]
    return expr_key


def _extract_loc_from_expr_key(expr_key):
    """Extract LOC part from expression key 'LOC:TCONS'."""
    if ':' in expr_key:
        return expr_key.split(':')[0]
    return expr_key


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

def write_module_outputs(module_id, module_name, expanded_isoforms,
                         ranking_df, rising_df, G, expr_filt,
                         term_matches, query_info,
                         gene_to_isoforms, gene_key_to_query, out_prefix):
    """Write all output files for one module.

    Parameters
    ----------
    module_id : str
        e.g. "M3"
    module_name : str
        Safe file-system name, e.g. "PPR_Proteins"
    expanded_isoforms : set
        Query IDs in the expanded module.
    ranking_df : pd.DataFrame
        Combined regulatory ranking (gene-level index: LOC IDs).
    rising_df : pd.DataFrame
        Rising hubs subset.
    G : nx.Graph
        Co-expression graph (gene-level nodes).
    expr_filt : pd.DataFrame
        Filtered expression matrix (gene-level).
    term_matches : dict
        query_id ? list of matched terms.
    query_info : dict
        query_id ? annotation info.
    gene_to_isoforms : dict
        gene_id (LOC*) ? list of isoform expression keys.
    gene_key_to_query : dict
        gene_id ? best query_id for annotation lookup.
    out_prefix : str
        Output file prefix.
    """
    safe_name = re.sub(r'[^a-zA-Z0-9_-]', '_', module_name)
    prefix = f"{out_prefix}_module_{module_id}_{safe_name}"

    # Build reverse mapping: gene ? query_ids and gene ? representative isoform info
    gene_to_queries = defaultdict(list)
    for gid, iso_keys in gene_to_isoforms.items():
        for ik in iso_keys:
            qid = gene_key_to_query.get(gid, gid)
            gene_to_queries[gid].append(qid)
    if not gene_to_queries:
        # Fallback: use gene IDs directly
        for gid in set(ranking_df.index):
            gene_to_queries[gid] = [gene_key_to_query.get(gid, gid)]

    # 1. Expanded nodes table
    nodes_path = f"{prefix}_expanded_nodes.tsv"
    node_rows = []
    for qid in sorted(expanded_isoforms):
        info = query_info.get(qid, {})
        matches = term_matches.get(qid, [])
        # Find which gene this belongs to in expression data
        gene_id = ''
        for gid, iso_keys in gene_to_isoforms.items():
            for ik in iso_keys:
                if qid in ik or (':' in ik and qid == ik.split(':')[1]):
                    gene_id = gid
                    break
            if gene_id:
                break
        if not gene_id:
            # Fallback: use LOC part if available
            gene_id = qid.split(':')[0] if ':' in qid else qid
        # Collect isoforms for this gene
        isoforms_of_gene = gene_to_isoforms.get(gene_id, [])
        node_rows.append({
            'query_id': qid,
            'gene_id': gene_id,
            'n_isoforms_in_expr': len(isoforms_of_gene),
            'isoforms_in_expr': ', '.join(isoforms_of_gene[:10]),
            'annotation': urllib.parse.unquote(
                info.get('desc', '')) if isinstance(info.get('desc'), str)
                else info.get('desc', ''),
            'pfam': ', '.join(info.get('pfam', [])),
            'go_ids': ', '.join(info.get('go', [])),
            'matched_terms': '; '.join(matches),
            'n_matched_terms': len(matches),
        })
    pd.DataFrame(node_rows).to_csv(nodes_path, sep='\t', index=False)
    print(f"  ? {nodes_path}  ({len(node_rows)} expanded query IDs)")

    # 2. Co-expression edges table (gene-level)
    if G.number_of_edges() > 0:
        edges_path = f"{prefix}_coexpression_edges.tsv"
        edge_rows = []
        for u, v, d in G.edges(data=True):
            # Annotate with key isoforms for each gene
            iso_u = gene_key_to_query.get(u, u)
            iso_v = gene_key_to_query.get(v, v)
            edge_rows.append({
                'source_gene': u,
                'target_gene': v,
                'source_best_query': iso_u,
                'target_best_query': iso_v,
                'rho': d.get('weight', 0),
                'pvalue': d.get('pvalue', 1),
            })
        pd.DataFrame(edge_rows).to_csv(edges_path, sep='\t', index=False)
        print(f"  ? {edges_path}  ({len(edge_rows)} edges)")

    # 3. Regulatory ranking table (gene-level, annotated with isoforms)
    rank_path = f"{prefix}_regulatory_ranking.tsv"
    cols = [
        'combined_score', 'rank_combined', 'rank_gwas', 'rank_change',
        'gwas_pvalue', 'neg_log10_p', 'evidence', 'has_gwas_signal',
        'degree_centrality', 'betweenness_centrality',
        'eigenvector_centrality', 'kME', 'clustering_coefficient',
        'z_gwas', 'z_betweenness', 'z_kME',
    ]
    avail_cols = [c for c in cols if c in ranking_df.columns]
    # Add gene annotation columns
    ranking_out = ranking_df[avail_cols].copy()
    ranking_out['query_ids'] = [', '.join(gene_to_queries.get(g, [g]))
                                 for g in ranking_df.index]
    ranking_out['n_isoforms'] = [len(gene_to_isoforms.get(g, []))
                                  for g in ranking_df.index]
    ranking_out['isoforms'] = [', '.join(gene_to_isoforms.get(g, [g])[:5])
                                for g in ranking_df.index]
    ranking_out.to_csv(rank_path, sep='\t', index=True, index_label='gene_id')
    print(f"  ? {rank_path}  ({len(ranking_out)} ranked genes)")
    if 'evidence' in ranking_out.columns:
        n_gwas = (ranking_out['evidence'] == 'gwas+expr').sum()
        n_expr = (ranking_out['evidence'] == 'expr_only').sum()
        print(f"     evidence: {n_gwas} gwas+expr  |  {n_expr} expr-only")

    # 4. Rising hubs (module-specific)
    if not rising_df.empty:
        rise_path = f"{prefix}_rising_hubs.tsv"
        rising_out = rising_df[avail_cols].copy()
        rising_out['query_ids'] = [', '.join(gene_to_queries.get(g, [g]))
                                    for g in rising_df.index]
        rising_out['n_isoforms'] = [len(gene_to_isoforms.get(g, []))
                                     for g in rising_df.index]
        rising_out['isoforms'] = [', '.join(gene_to_isoforms.get(g, [g])[:5])
                                   for g in rising_df.index]
        rising_out.to_csv(rise_path, sep='\t', index=True, index_label='gene_id')
        print(f"  ? {rise_path}  ({len(rising_out)} rising hubs)")
        # Print top rising hubs
        print(f"\n  {'='*70}")
        print(f"  Module {module_id} ({module_name}) : Top Rising Hubs")
        print(f"  {'='*70}")
        for idx, (gene, row) in enumerate(rising_out.head(15).iterrows()):
            pval = row.get('gwas_pvalue', 1.0)
            combined = row.get('combined_score', 0)
            kme = row.get('kME', 0)
            btwn = row.get('betweenness_centrality', 0)
            qids = row.get('query_ids', gene)
            print(f"  {idx+1:2d}. {gene:30s}  "
                  f"GWAS p={pval:.2e}  ?  combined={combined:.3f}  "
                  f"(?rank=+{int(row['rank_change'])})  "
                  f"kME={kme:.3f}  btwn={btwn:.3f}  "
                  f"isoforms={qids[:60]}")


def write_cross_module_rising_hubs(all_rising, out_prefix):
    """Write a consolidated cross-module rising hub report.

    Parameters
    ----------
    all_rising : list of (module_id, module_name, pd.DataFrame)
    out_prefix : str
    """
    if not all_rising:
        return

    all_rows = []
    for mod_id, mod_name, rising_df in all_rising:
        for iso, row in rising_df.iterrows():
            all_rows.append({
                'module_id': mod_id,
                'module_name': mod_name,
                'isoform': iso,
                'combined_score': row.get('combined_score', 0),
                'rank_combined': row.get('rank_combined', 0),
                'rank_gwas': row.get('rank_gwas', 0),
                'rank_change': int(row.get('rank_change', 0)),
                'gwas_pvalue': row.get('gwas_pvalue', 1.0),
                'kME': row.get('kME', 0),
                'betweenness_centrality': row.get('betweenness_centrality', 0),
            })

    df = pd.DataFrame(all_rows)
    df = df.sort_values('rank_change', ascending=False)

    path = f"{out_prefix}_rising_hubs.tsv"
    df.to_csv(path, sep='\t', index=False)
    print(f"\n{'='*70}")
    print(f"  CROSS-MODULE RISING HUBS ? {path}  ({len(df)} hubs)")
    print(f"{'='*70}")
    for idx, row in df.head(30).iterrows():
        print(f"  {row['module_id']:4s}  {row['isoform']:45s}  "
              f"GWAS p={row['gwas_pvalue']:.2e}  ?  "
              f"combined={row['combined_score']:.3f}  "
              f"(?rank=+{row['rank_change']})  "
              f"kME={row['kME']:.3f}")


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="ISONET_LAYER2 — expression-based network expansion and regulatory ranking",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Input files
    parser.add_argument('--modules-tsv', required=True,
                        help='Co-function modules TSV from ISONET_LAYER1')
    parser.add_argument('--nodes-tsv', required=True,
                        help='Co-function nodes TSV from ISONET_LAYER1')
    parser.add_argument('--edges-tsv', required=True,
                        help='Co-function edges TSV from ISONET_LAYER1')
    parser.add_argument('--annotations', required=True,
                        help='Full eggNOG-mapper annotations file')
    parser.add_argument('--expression', required=True,
                        help='TPM expression matrix (Pop_isoform_TPM_by_compareGTF.txt)')
    parser.add_argument('--gtf', default=None,
                        help='StringTie GTF for TCONS?LOC mapping')
    parser.add_argument('--qqman', required=True,
                        help='QQMan GWAS summary file for p-value lookup')
    parser.add_argument('--go-tab', dest='go_tab', default=None,
                        help='go.tab for GO ID ? description mapping')

    # Output
    parser.add_argument('--out-prefix', dest='out_prefix',
                        default='layer2_expression',
                        help='Prefix for output files')

    # Module filtering : by ID, size, or AI-assigned confidence level
    parser.add_argument('--min-module-size', dest='min_module_size',
                        type=int, default=3,
                        help='Minimum module size to expand and analyze')
    parser.add_argument('--modules', nargs='+', default=None,
                        help='Specific module IDs to analyze (e.g., M3 M9). '
                             'If omitted, filtered by --ai-confidence and --min-module-size.')
    parser.add_argument('--ai-analysis-json', dest='ai_analysis_json', default=None,
                        help='Path to _network_ai_analysis.json from ISONET_LAYER1 '
                             '(contains per-module AI confidence scores: high/medium/low). '
                             'Used with --ai-confidence to select only high-confidence modules.')
    parser.add_argument('--ai-confidence', dest='ai_confidence',
                        default='high',
                        choices=['high', 'high,medium', 'high,medium,low', 'all'],
                        help='Minimum AI confidence level for module selection. '
                             '"high" = only high-confidence modules (default : most reliable). '
                             '"high,medium" = high + medium. '
                             '"high,medium,low" = all modules with any AI assessment. '
                             '"all" = skip AI filtering entirely. '
                             'Ignored if --modules is explicitly provided.')

    # Expansion parameters
    parser.add_argument('--min-shared-terms', dest='min_shared_terms',
                        type=int, default=2,
                        help='Minimum shared functional terms for annotation expansion')
    parser.add_argument('--max-expansion', dest='max_expansion',
                        type=int, default=500,
                        help='Maximum expanded module size')

    # Multi-hop expansion: follow co-expression links from original module
    # members to capture functionally linked genes missed by annotation alone.
    parser.add_argument('--expansion-hops', dest='expansion_hops',
                        type=int, default=1, choices=[1, 2, 3],
                        help='Expansion depth: 1=annotation only, '
                             '2=also include co-expressed neighbors of module '
                             'members, 3=include neighbors of neighbors. '
                             'Higher hops capture more regulatory context but '
                             'increase runtime.')
    parser.add_argument('--expansion-link-rho', dest='expansion_link_rho',
                        type=float, default=0.6,
                        help='Minimum |Spearman ?| for a co-expression edge '
                             'to trigger inclusion during multi-hop expansion. '
                             'Lower values include more genes (broader network); '
                             'higher values are more stringent.')

    # TF-target expansion: for transcription factor modules, correlate TF seed
    # expression against the entire TPM matrix to find candidate target genes.
    parser.add_argument('--tf-min-corr', dest='tf_min_corr',
                        type=float, default=0.3,
                        help='Minimum |Pearson r| for a gene to be considered '
                             'a TF target in expression correlation screen. '
                             'Lower values include more candidate targets.')
    parser.add_argument('--tf-top-n', dest='tf_top_n',
                        type=int, default=300,
                        help='Maximum target genes to return per TF seed.')

    # Co-expression parameters
    parser.add_argument('--corr-threshold', dest='corr_threshold',
                        type=float, default=0.7,
                        help='Minimum |Spearman ?| for co-expression edge')
    parser.add_argument('--fdr-threshold', dest='fdr_threshold',
                        type=float, default=0.05,
                        help='FDR threshold for co-expression edges')
    parser.add_argument('--min-tpm', dest='min_tpm',
                        type=float, default=0.1,
                        help='Minimum mean TPM to retain an isoform')

    # Scoring weights
    parser.add_argument('--gwas-weight', dest='gwas_weight',
                        type=float, default=0.3,
                        help='Weight for GWAS evidence in combined score')
    parser.add_argument('--expr-weight', dest='expr_weight',
                        type=float, default=0.7,
                        help='Weight for expression evidence in combined score')
    parser.add_argument('--expression-only', dest='expression_only',
                        action='store_true', default=False,
                        help='Score genes purely by expression centrality '
                             '(GWAS weight = 0).  Use to discover regulatory '
                             'hubs that GWAS cannot detect : common for '
                             'essential genes under purifying selection, '
                             'post-transcriptional regulators, and redundant '
                             'pathway members.')

    # Rising hub detection
    parser.add_argument('--label-density', dest='label_density',
                        type=str, default='auto',
                        choices=['sparse', 'normal', 'dense', 'all', 'auto'],
                        help='Label density for the R network plot. '
                             '"sparse" = GWAS seeds + highlights only. '
                             '"normal" = + top 10 rising hubs + top 5 by degree. '
                             '"dense" = + neighbors + top 20 rising + kME > 0.7. '
                             '"all" = every node labeled (crowded). '
                             '"auto" = dense for < 100 nodes, normal otherwise. '
                             'Set via LABEL_DENSITY env var for the R script.')
    parser.add_argument('--force-genes', dest='force_genes',
                        type=str, default='',
                        help='Comma-separated gene IDs (LOC*/TCONS*) to force-include '
                             'in every expanded module. These genes bypass annotation '
                             'scoring and are always loaded for co-expression. '
                             'Use for genes of interest that share module terms but '
                             'fall below the top-N cutoff.')
    parser.add_argument('--min-rank-rise', dest='min_rank_rise',
                        type=int, default=10,
                        help='Minimum rank improvement for rising hub status')
    parser.add_argument('--min-gwas-pvalue', dest='min_gwas_pvalue',
                        type=float, default=0.01,
                        help='Max GWAS p-value to be a "rising hub"')

    args = parser.parse_args()

    t_total = time.time()

    # ------------------------------------------------------------------
    # 1. Load co-function modules and nodes
    # ------------------------------------------------------------------
    print(f"[1/7] Loading co-function modules ?")
    modules_df = pd.read_csv(args.modules_tsv, sep='\t')
    nodes_df = pd.read_csv(args.nodes_tsv, sep='\t')

    # Build nodes dict: isoform_id ? annotation fields
    nodes = {}
    for _, row in nodes_df.iterrows():
        nodes[row['Isoform']] = row.to_dict()

    # Parse module isoform lists
    modules = {}
    # Load AI analysis first (needed for module names + confidence)
    ai_module_conf = {}   # module_id ? confidence
    ai_module_names = {}  # module_id ? human-readable AI name
    if args.ai_analysis_json and os.path.exists(args.ai_analysis_json):
        with open(args.ai_analysis_json) as f:
            ai_data = json.load(f)
        for m in ai_data.get('modules', []):
            mid = m['module_id']
            ai_module_conf[mid] = m.get('confidence', 'low')
            ai_name = m.get('name', '')
            if ai_name:
                ai_module_names[mid] = re.sub(r'[^a-zA-Z0-9_-]', '_', ai_name)[:40]
        n_high = sum(1 for v in ai_module_conf.values() if v == 'high')
        n_med  = sum(1 for v in ai_module_conf.values() if v == 'medium')
        n_low  = sum(1 for v in ai_module_conf.values() if v == 'low')
        print(f"  AI confidence: {n_high} high, {n_med} medium, {n_low} low")

    module_names = {}
    for _, row in modules_df.iterrows():
        mid = row['module_id']
        iso_list = [i.strip() for i in str(row['isoforms']).split(',') if i.strip()]
        modules[mid] = iso_list
        # Use AI name if available, otherwise fall back to top PFAM term
        if mid in ai_module_names:
            module_names[mid] = ai_module_names[mid]
        else:
            top_terms = str(row.get('top_terms', ''))
            if top_terms and top_terms != 'nan':
                # Extract best term: prefer PFAM > GO > COG
                best = ''
                for term in top_terms.split(';'):
                    term = term.strip()
                    if term.startswith('PFAM:'):
                        best = term.split(':')[1].strip()
                        break
                    elif term.startswith('GO_description:') and not best:
                        best = term.split(':')[1].strip()
                if not best:
                    first = top_terms.split(';')[0].split(':')[-1].strip().split()[0]
                    best = first[:30]
                module_names[mid] = re.sub(r'[^a-zA-Z0-9_-]', '_', best)[:30]
            else:
                module_names[mid] = mid

    print(f"  {len(modules)} modules loaded, {len(nodes)} nodes")

    # Determine which confidence levels to include
    if args.ai_confidence == 'all':
        allowed_conf = {'high', 'medium', 'low'}
    else:
        allowed_conf = set(args.ai_confidence.split(','))

    # Filter modules
    if args.modules:
        # Explicit module list : bypass AI filtering
        target_modules = {m for m in args.modules if m in modules}
    else:
        # Auto-select: filter by size AND AI confidence
        target_modules = set()
        for m, isos in modules.items():
            if len(isos) < args.min_module_size:
                continue
            # If AI confidence data is available, filter by it
            if ai_module_conf:
                conf = ai_module_conf.get(m, 'low')
                if conf not in allowed_conf:
                    continue
            target_modules.add(m)

    target_modules = sorted(target_modules,
                           key=lambda m: len(modules[m]), reverse=True)

    # Build summary with confidence labels
    print(f"  Analyzing {len(target_modules)} modules "
          f"(ai_confidence ? {args.ai_confidence}, "
          f"min_size ? {args.min_module_size}):")
    for m in target_modules:
        conf = ai_module_conf.get(m, '?')
        print(f"    {m}: {len(modules[m])} isoforms  [{conf}]")

    # ------------------------------------------------------------------
    # 2. Build full annotation index
    # ------------------------------------------------------------------
    print(f"\n[2/7] Building annotation index ?")
    ann_idx, query_info = build_annotation_index(args.annotations)

    # Load GO map if available
    go_map = {}
    if args.go_tab and os.path.exists(args.go_tab):
        print(f"[GO] Loading GO descriptions: {args.go_tab}")
        with open(args.go_tab, 'r') as f:
            _ = f.readline()  # skip header
            for line in f:
                parts = line.rstrip('\n').split('\t')
                if len(parts) >= 2:
                    go_map[parts[0]] = parts[1]
        print(f"[GO]   {len(go_map):,} GO terms loaded")

    # ------------------------------------------------------------------
    # 3. Load GWAS p-values
    # ------------------------------------------------------------------
    print(f"\n[3/7] Loading GWAS p-values ?")
    gwas_pvalues = load_gwas_pvalues(args.qqman)

    # ------------------------------------------------------------------
    # 4. Parse GTF for TCONS?LOC mapping (optional)
    # ------------------------------------------------------------------
    gtf_tcons_to_loc = {}
    if args.gtf and os.path.exists(args.gtf):
        print(f"\n[4/7] Parsing GTF for ID mapping ?")
        t0 = time.time()
        with open(args.gtf, 'r') as fh:
            for line in fh:
                if line.startswith('#') or '\ttranscript\t' not in line:
                    continue
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 9:
                    continue
                # Parse attributes
                attrs = {}
                for m in re.finditer(r'(\S+)\s+"([^"]*)"', parts[8]):
                    attrs[m.group(1)] = m.group(2)
                tid = attrs.get('transcript_id')
                gene_name = attrs.get('gene_name', '')
                if tid and gene_name:
                    gtf_tcons_to_loc[tid] = gene_name
        print(f"  [GTF] {len(gtf_tcons_to_loc):,} TCONS?LOC mappings "
              f"({time.time() - t0:.1f}s)")
    else:
        print(f"\n[4/7] Skipping GTF (not provided)")

    # ------------------------------------------------------------------
    # 5. Pre-scan expression file for available gene/isoform IDs
    #    (used as a bonus/tiebreaker in module expansion so
    #     expression-relevant genes aren't crowded out)
    # ------------------------------------------------------------------
    print(f"\n[5/7] Scanning expression data for available IDs ?")
    t_scan = time.time()
    expr_available_ids = _scan_expression_ids(args.expression)
    print(f"  {len(expr_available_ids):,} unique IDs found "
          f"({time.time() - t_scan:.1f}s)")

    # ------------------------------------------------------------------
    # 6. For each module: extract terms, expand, build co-expression
    # ------------------------------------------------------------------
    print(f"\n[6/7] Processing modules ?")

    all_rising = []    # (module_id, module_name, rising_df)

    for mid in target_modules:
        mod_isos = modules[mid]
        mod_name = module_names.get(mid, mid)

        # Keep only isoforms that are in nodes dict
        mod_isos_in_nodes = [i for i in mod_isos if i in nodes]
        if len(mod_isos_in_nodes) < 2:
            print(f"\n  [{mid}] SKIP: <2 isoforms with node data")
            continue

        print(f"\n  {'?'*60}")
        print(f"  [{mid}] {mod_name} : {len(mod_isos_in_nodes)} isoforms")
        print(f"  {'?'*60}")

        # 5a. Extract characterizing terms
        mod_terms = extract_module_terms(mod_isos_in_nodes, nodes)
        n_terms = sum(len(v) for v in mod_terms.values())
        print(f"  [TERMS] {n_terms} characterizing terms extracted:")
        for cat, terms in mod_terms.items():
            if terms:
                top = terms[:5]
                print(f"    {cat}: {', '.join(t for t, _ in top)}")

        # 5b. Expansion: annotation-based vs TF-target correlation.
        #
        # For **TF modules** the annotation expansion is actively harmful:
        # searching on "Myb_DNA-binding" pulls in every MYB in the genome
        # regardless of whether they regulate chlorophyll.  Instead, we
        # correlate each TF seed's expression vector against the entire TPM
        # matrix to find its downstream targets — genes whose transcript
        # levels track with the TF in vivo.
        #
        # For **non-TF modules**, the original annotation expansion
        # (shared PFAM/GO/KEGG terms) is correct — it captures functionally
        # related proteins like PPR editing factors or kinase cascades.
        is_tf = _is_tf_module(mod_terms)
        expanded = set(mod_isos_in_nodes)  # always keep original module members
        term_matches = {}                  # populated only by annotation expansion

        # Force-include user-specified genes (e.g. --force-genes LOC111886547)
        if args.force_genes:
            forced = {g.strip() for g in args.force_genes.split(',') if g.strip()}
            expanded |= forced
            if forced:
                n_forced_new = len(forced - set(mod_isos_in_nodes))
                if n_forced_new > 0:
                    print(f"  [FORCE] {n_forced_new} gene(s) force-included: "
                          f"{', '.join(sorted(forced - set(mod_isos_in_nodes)))}")

        if is_tf:
            print(f"  [TF-TARGET] Module {mid} is a TF module, "
                  f"skipping annotation expansion (avoids pulling in "
                  f"all Myb/bHLH/etc. proteins).")
            if n_terms > 0:
                print(f"  [TF-TARGET]   Dropped {n_terms} annotation terms: "
                      f"{', '.join(t for cat, terms in mod_terms.items() for t, _ in terms[:3])}")

            # Separate seed genes into true-TF and non-TF seeds.
            # Only true TFs (with known TF PFAM domains) get target expansion;
            # non-TF seeds are simply carried through as GWAS-validated members.
            tf_seed_genes = []   # list of (gene_id, pfam_label) for display
            non_tf_seeds = set()
            for iso in mod_isos_in_nodes:
                nd = nodes.get(iso, {})
                iso_pfam = set(str(nd.get('PFAM', '')).split(','))
                iso_pfam = {t.strip() for t in iso_pfam if t.strip() and t.strip() != '-'}
                is_tf_seed = bool(iso_pfam & _TF_PFAM_TERMS)

                # Resolve to gene ID
                if iso.startswith('LOC'):
                    gene_id = iso
                elif iso.startswith('TCONS') and gtf_tcons_to_loc:
                    gene_id = gtf_tcons_to_loc.get(iso)
                else:
                    gene_id = None

                if is_tf_seed and gene_id:
                    # Avoid duplicates (TCONS isoform may map to same LOC)
                    existing = {g for g, _ in tf_seed_genes}
                    if gene_id not in existing:
                        pfam_label = ', '.join(sorted(iso_pfam & _TF_PFAM_TERMS))
                        tf_seed_genes.append((gene_id, pfam_label))
                        print(f"  [TF-TARGET]   TF seed: {iso} → {gene_id} "
                              f"(PFAM: {pfam_label})")
                elif not is_tf_seed:
                    non_tf_seeds.add(iso)
                    print(f"  [TF-TARGET]   non-TF seed: {iso} — kept as member, not used for target search")

            if not tf_seed_genes:
                print(f"  [TF-TARGET]   no valid TF seeds found — "
                      f"module members may all be co-associated non-TF genes")

            # --- Per-TF sub-module expansion -----------------------------------
            # Each TF seed gets its own expansion so its target network is
            # visually distinct from other TFs in the same module.
            for tf_idx, (tf_gene_id, tf_pfam_label) in enumerate(tf_seed_genes):
                # Unique sub-module ID for output files
                sub_mid = f"{mid}_{tf_gene_id}"
                # Clean label: parent module name + "targets of" + TF gene
                sub_name = f"{mod_name} targets of {tf_gene_id}"
                # Keep filename-safe
                sub_name = re.sub(r'[^a-zA-Z0-9_ -]', '', sub_name)[:80]

                print(f"\n  {'?'*60}")
                print(f"  [TF-SUB] {sub_mid} — {tf_gene_id} ({tf_pfam_label})")
                print(f"  {'?'*60}")

                # Build expanded set: this TF + non-TF co-associates + force genes
                sub_expanded = {tf_gene_id} | non_tf_seeds
                if args.force_genes:
                    forced = {g.strip() for g in args.force_genes.split(',') if g.strip()}
                    sub_expanded |= forced

                # Find targets for just this one TF
                sub_targets, _ = find_tf_targets(
                    args.expression, {tf_gene_id},
                    expr_available_ids=expr_available_ids,
                    min_corr=getattr(args, 'tf_min_corr', 0.3),
                    top_n=getattr(args, 'tf_top_n', 300),
                )
                if sub_targets:
                    sub_expanded |= sub_targets
                    print(f"  [TF-TARGET]   {len(sub_targets)} targets for {tf_gene_id} "
                          f"(now {len(sub_expanded)} total)")
                else:
                    print(f"  [TF-TARGET]   no targets found for {tf_gene_id}")

                if len(sub_expanded) < 3:  # need TF + at least 2 targets
                    print(f"  [TF-SUB] SKIP: insufficient expansion for {tf_gene_id}")
                    continue

                # ── Run shared downstream for this TF sub-module ──────────────
                _process_one_expansion(
                    sub_mid, sub_name, sub_expanded,
                    mod_isos_in_nodes,
                    args, expr_available_ids, gtf_tcons_to_loc,
                    gwas_pvalues, ann_idx, query_info, go_map,
                    all_rising,
                )

            # TF modules handled — skip the shared single-expansion path below
            continue

        else:
            # Non-TF module: standard annotation-based expansion
            if n_terms == 0:
                print(f"  [{mid}] SKIP: no specific terms to expand with")
                continue

            new_from_ann, term_matches = expand_module(
                mod_terms, ann_idx, query_info, go_map=go_map,
                min_shared_terms=args.min_shared_terms,
                max_expansion=args.max_expansion,
                expr_available_ids=expr_available_ids,
            )
            expanded |= new_from_ann
            n_new = len(expanded - set(mod_isos_in_nodes))
            print(f"  [EXPAND] {len(expanded)} isoforms matched "
                  f"({n_new} new, beyond original {len(mod_isos_in_nodes)})")

        _process_one_expansion(
            mid, mod_name, expanded, mod_isos_in_nodes,
            args, expr_available_ids, gtf_tcons_to_loc,
            gwas_pvalues, term_matches, query_info, go_map,
            all_rising,
        )

    # ------------------------------------------------------------------
    # 7. Cross-module rising hub report
    # ------------------------------------------------------------------
    print(f"\n[7/7] Cross-module rising hub report ?")
    if all_rising:
        write_cross_module_rising_hubs(all_rising, args.out_prefix)
    else:
        print("  No rising hubs detected across modules.")

    # ------------------------------------------------------------------
    # 8. Summary
    # ------------------------------------------------------------------
    t_elapsed = time.time() - t_total
    print(f"\n[8/8] Done in {t_elapsed:.1f}s total")
    print(f"  Output prefix: {args.out_prefix}")


if __name__ == '__main__':
    main()
