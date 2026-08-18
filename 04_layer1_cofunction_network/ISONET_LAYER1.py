#!/usr/bin/env python3
"""
ISONET_LAYER1.py — Multi-metabolite GWAS to annotated isoform table
                            (Layer 1: GWAS → annotation → co-function network)

Batch-processes QQMan GWAS summary files (one per metabolite), extracts
significant isoforms (filtered by p-value), maps each isoform to its parent
gene via StringTie GTF, annotates via eggNOG-mapper, resolves GO term
descriptions, optionally enriches results with AI, and builds a co-function
network of functionally related isoforms.

Dual-provider AI support (OpenAI + DeepSeek):
    - --thinking flag auto-detects provider: reasoning for OpenAI, thinking for DeepSeek
    - --reasoning-effort flag (low / medium / high, default: medium)
    - Fully backward-compatible: DeepSeek models still use thinking {type: enabled}

Output:
    - {prefix}_significant_isoforms.tsv : consolidated annotated table
    - {prefix}_summary_stats.txt        : per-metabolite match statistics

Columns (output table):
    Metabolites   — metabolite/phenotype name (from filename)
    Chr           — chromosome in Lsat_Salinas_v11_chrN format
    Gene_start    — parent gene start coordinate (from GFF)
    Gene_end      — parent gene end coordinate (from GFF)
    ref_Gene      — parent reference gene (LOC ID)
    Isoform       — significant transcript/isoform (TCONS or LOC ID)
    Pvalue        — GWAS p-value
    Annotation    — eggNOG Description (decoded)
    GO_description— top resolved GO term descriptions (human-readable)
    PFAM          — PFAM domain names
    AI_link       — (optional) AI-generated hypothesis linking isoform to metabolite
    AI_confidence — (optional) AI-assigned confidence: high | medium | low | none

Usage:
    python3 ISONET_LAYER1.py \\
        --qqman-dir /path/to/metabolites/groupB/ \\
        --qqman-pattern "*.pheno.qqman" \\
        --gtf pop240GTF2ref.combined.gtf \\
        --gff Lsat_Salinas_v11_genomic.gff \\
        --annotations isoform_and_gene.emapper.annotations \\
        --go-tab go.tab \\
        --pvalue-threshold 1e-5 \\
        --out-prefix metabolites_significant

AI enrichment (OpenAI with reasoning):
    python3 ISONET_LAYER1.py ... \\
        --api-key YOUR_OPENAI_KEY \\
        --model gpt-5.5 \\
        --api-base https://api.openai.com/v1 \\
        --thinking \\
        --reasoning-effort medium \\
        --ai-system "You are a plant metabolomics expert. ..."

AI enrichment (DeepSeek, backward-compatible):
    python3 ISONET_LAYER1.py ... \\
        --api-key YOUR_DEEPSEEK_KEY \\
        --model deepseek-v4-pro \\
        --api-base https://api.deepseek.com/v1 \\
        --thinking
"""

import argparse
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request
from collections import OrderedDict
from math import log10

import pandas as pd


# ---------------------------------------------------------------------------
# Regex for extracting gene identifiers from QQMan entry names
# ---------------------------------------------------------------------------
RE_LOC = re.compile(r'LOC\d+')
RE_TCONS = re.compile(r'TCONS_\d+')


def extract_gene_id(raw_id):
    """Extract a canonical gene/transcript identifier from a QQMan row ID.

    Parameters
    ----------
    raw_id : str
        Raw entry from QQMan column 1, e.g. ``1_43133:LOC111905264``,
        ``TCONS_00004469``, or long-form with SV tags.

    Returns
    -------
    str
        The LOC or TCONS identifier, or the original string if neither found.
    """
    m = RE_LOC.search(raw_id)
    if m:
        return m.group(0)
    m = RE_TCONS.search(raw_id)
    if m:
        return m.group(0)
    return raw_id


# ---------------------------------------------------------------------------
# GTF parsing: TCONS → parent gene + chromosome
# ---------------------------------------------------------------------------

def parse_gtf_transcripts(gtf_path):
    """Parse StringTie combined GTF, build TCONS→{chr,start,end,gene_name} map.

    Only ``transcript`` feature lines are processed.

    Parameters
    ----------
    gtf_path : str
        Path to combined GTF (e.g. pop240GTF2ref.combined.gtf).

    Returns
    -------
    dict
        Keyed by transcript_id (TCONS_*) → dict with keys:
        ``chr``, ``start`` (int), ``end`` (int), ``gene_name`` (LOC* or other).
    """
    print(f"[GTF] Parsing transcript entries from: {gtf_path}")
    tcons_map = OrderedDict()
    attr_re = re.compile(r'(\S+)\s+"([^"]*)"')

    opener = urllib.request if gtf_path.endswith('.gz') else open
    open_kwargs = {}
    if gtf_path.endswith('.gz'):
        import gzip
        fh = gzip.open(gtf_path, 'rt')
    else:
        fh = open(gtf_path, 'r')

    try:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            if parts[2] != 'transcript':
                continue

            # Parse attributes
            attrs = {}
            for m in attr_re.finditer(parts[8]):
                attrs[m.group(1)] = m.group(2)

            tid = attrs.get('transcript_id')
            if tid is None:
                continue
            gene_name = attrs.get('gene_name', '')
            tcons_map[tid] = {
                'chr': parts[0],
                'start': int(parts[3]),
                'end': int(parts[4]),
                'gene_name': gene_name,
            }
    finally:
        fh.close()

    n_with_gene = sum(1 for v in tcons_map.values() if v['gene_name'])
    print(f"[GTF]   {len(tcons_map):,} transcripts loaded"
          f"  ({n_with_gene:,} with gene_name)")
    return tcons_map


# ---------------------------------------------------------------------------
# GFF parsing: gene entries + chromosome name mapping
# ---------------------------------------------------------------------------

def _parse_gff_attributes(attr_string):
    """Parse GFF3 column-9 attribute string into a dict."""
    attrs = {}
    for part in attr_string.split(';'):
        part = part.strip()
        if '=' in part:
            key, val = part.split('=', 1)
            attrs[key] = val
    return attrs


def parse_gff(gff_path):
    """Parse RefSeq GFF3 — returns gene dict and chromosome name map.

    Parameters
    ----------
    gff_path : str
        Path to RefSeq GFF3 file (e.g. Lsat_Salinas_v11_genomic.gff).

    Returns
    -------
    gene_map : dict
        Keyed by LOC ID → dict with keys:
        ``chr`` (Lsat_Salinas_v11_chrN), ``start`` (int), ``end`` (int),
        ``description``, ``biotype``.
    chr_map : dict
        NC_* → Lsat_Salinas_v11_chrN mapping (from ``region`` lines).
    """
    print(f"[GFF] Parsing gene and region entries from: {gff_path}")
    gene_map = OrderedDict()
    chr_map = {}
    region_count = 0

    import gzip as _gzip
    opener = _gzip.open if gff_path.endswith('.gz') else open

    with opener(gff_path, 'rt') if gff_path.endswith('.gz') else open(gff_path, 'r') as fh:
        for line in fh:
            if line.startswith('#') or line.strip() == '':
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            feat_type = parts[2]
            chr_id = parts[0]

            if feat_type == 'region':
                attrs = _parse_gff_attributes(parts[8])
                name = attrs.get('Name', attrs.get('chromosome', ''))
                if name and name not in ('Unknown', ''):
                    try:
                        chr_num = int(name)
                        chr_map[chr_id] = f'Lsat_Salinas_v11_chr{chr_num}'
                    except ValueError:
                        chr_map[chr_id] = f'Lsat_Salinas_v11_chr{name}'
                region_count += 1

            elif feat_type == 'gene':
                attrs = _parse_gff_attributes(parts[8])
                gene_id = attrs.get('ID', '').replace('gene-', '')
                # Ensure it's a LOC identifier
                if not gene_id.startswith('LOC'):
                    gene_id = attrs.get('Name', gene_id)
                if not gene_id.startswith('LOC'):
                    continue  # skip non-LOC gene entries

                gene_map[gene_id] = {
                    'chr': chr_id,       # NC_* — resolved later via chr_map
                    'start': int(parts[3]),
                    'end': int(parts[4]),
                    'description': attrs.get('description', attrs.get('product', '')),
                    'biotype': attrs.get('gene_biotype', ''),
                }

    # Resolve chromosome names in gene_map
    resolved = 0
    for gid, rec in gene_map.items():
        nc = rec['chr']
        if nc in chr_map:
            rec['chr'] = chr_map[nc]
            resolved += 1
        elif nc.startswith('NW_') or nc.startswith('NT_'):
            rec['chr'] = nc  # keep unplaced scaffold name as-is
        else:
            # Try numeric extraction from NC_*
            m = re.match(r'NC_\d+\.\d+', nc)
            if m:
                rec['chr'] = nc  # couldn't map, keep as-is

    print(f"[GFF]   {len(gene_map):,} LOC gene entries  |  "
          f"{region_count} regions  |  {resolved}/{len(gene_map)} chr names resolved")
    return gene_map, chr_map


# ---------------------------------------------------------------------------
# eggNOG-mapper annotation loading
# ---------------------------------------------------------------------------

def load_emapper(emapper_path):
    """Parse eggNOG-mapper annotations file into a dict keyed by query ID.

    Parameters
    ----------
    emapper_path : str
        Path to eggNOG-mapper annotations file.

    Returns
    -------
    dict
        Keyed by query ID (TCONS_* or LOC*) → dict of annotation fields.
    """
    print(f"[EMAPPER] Loading annotations: {emapper_path}")
    annot = OrderedDict()

    with open(emapper_path, 'r') as fh:
        header_line = None
        for line in fh:
            if line.startswith('##'):
                continue
            header_line = line.strip()
            break
        if header_line is None:
            raise ValueError(f"Could not find header line in: {emapper_path}")

        cols = header_line.lstrip('#').split('\t')
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 8:
                continue
            qid = parts[0]
            record = {}
            for j, col_name in enumerate(cols):
                record[col_name] = parts[j] if j < len(parts) else ''
            annot[qid] = record

    n_tcons = sum(1 for k in annot if k.startswith('TCONS'))
    n_loc = sum(1 for k in annot if k.startswith('LOC'))
    print(f"[EMAPPER]   {len(annot):,} records  (TCONS: {n_tcons:,}  |  LOC: {n_loc:,})")
    return annot


# ---------------------------------------------------------------------------
# GO description loading
# ---------------------------------------------------------------------------

def load_go_map(go_tab_path):
    """Load go.tab (GO ID → description) from eggNOG-mapper data.

    Parameters
    ----------
    go_tab_path : str
        Path to go.tab (tab-separated: GO, Description, level).

    Returns
    -------
    dict
        Keyed by GO ID (e.g. ``GO:0003700``) → description string.
    """
    print(f"[GO] Loading GO descriptions: {go_tab_path}")
    go_map = {}
    with open(go_tab_path, 'r') as fh:
        _ = fh.readline()  # skip header
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 2:
                go_map[parts[0]] = parts[1]
    print(f"[GO]   {len(go_map):,} GO term descriptions loaded")
    return go_map


# ---------------------------------------------------------------------------
# QQMan processing
# ---------------------------------------------------------------------------

def process_qqman(qqman_path, metabolite_name, pvalue_threshold,
                  tcons_map, gff_gene_map, emapper, go_map):
    """Process one metabolite QQMan file → list of annotated hit rows.

    Parameters
    ----------
    qqman_path : str
        Path to a single QQMan file (tab-sep, cols: ID CHR POS PVALUE).
    metabolite_name : str
        Metabolite/phenotype name (derived from filename).
    pvalue_threshold : float
        Maximum p-value for inclusion.
    tcons_map : dict
        GTF TCONS → {chr, start, end, gene_name} mapping.
    gff_gene_map : dict
        GFF LOC → {chr, start, end, description, biotype} mapping.
    emapper : dict
        eggNOG-mapper annotations keyed by query ID.
    go_map : dict
        GO ID → human-readable description mapping.

    Returns
    -------
    list of dict
        Each dict has keys: Metabolites, Chr, Gene_start, Gene_end,
        ref_Gene, Isoform, Pvalue, Annotation, GO_description, PFAM.
    """
    # Read QQMan
    df = pd.read_csv(
        qqman_path, sep='\t', header=None,
        names=['original_id', 'chr', 'pos', 'pvalue'],
        dtype={'original_id': str, 'chr': str, 'pos': int, 'pvalue': float}
    )

    # Extract canonical gene/isoform ID
    df['gene_id'] = df['original_id'].apply(extract_gene_id)

    # Deduplicate: keep lowest p-value per gene_id
    df = df.loc[df.groupby('gene_id')['pvalue'].idxmin()].copy()

    # Filter by p-value threshold
    df = df[df['pvalue'] <= pvalue_threshold].copy()

    if len(df) == 0:
        return []

    # Sort by p-value
    df = df.sort_values('pvalue', ascending=True).reset_index(drop=True)

    rows = []
    n_tcons_found = 0
    n_loc_found = 0
    n_annot_found = 0

    for _, row in df.iterrows():
        gid = row['gene_id']
        pval = row['pvalue']

        # --- Resolve coordinates and parent gene ---
        if gid.startswith('TCONS'):
            isoform_id = gid
            # Look up TCONS in GTF for chr and parent gene
            tcons_info = tcons_map.get(gid)
            if tcons_info:
                chr_name = tcons_info['chr']  # Lsat_Salinas_v11_chrN (already resolved)
                ref_gene = tcons_info['gene_name']
                n_tcons_found += 1

                # Get gene coordinates from GFF (parent LOC gene)
                if ref_gene and ref_gene in gff_gene_map:
                    ginfo = gff_gene_map[ref_gene]
                    gene_start = ginfo['start']
                    gene_end = ginfo['end']
                else:
                    # Fallback to transcript coordinates from GTF
                    gene_start = tcons_info['start']
                    gene_end = tcons_info['end']
            else:
                # TCONS not in GTF — use QQMan data
                chr_name = _chr_num_to_lsat(row['chr'])
                ref_gene = ''
                gene_start = row['pos']
                gene_end = row['pos']

        elif gid.startswith('LOC'):
            isoform_id = gid
            ref_gene = gid
            # Get coordinates from GFF
            if gid in gff_gene_map:
                ginfo = gff_gene_map[gid]
                chr_name = ginfo['chr']
                gene_start = ginfo['start']
                gene_end = ginfo['end']
                n_loc_found += 1
            else:
                chr_name = _chr_num_to_lsat(row['chr'])
                gene_start = row['pos']
                gene_end = row['pos']
        else:
            # Unrecognized ID format
            isoform_id = gid
            ref_gene = ''
            chr_name = _chr_num_to_lsat(row['chr'])
            gene_start = row['pos']
            gene_end = row['pos']

        # --- Look up annotations ---
        annotation = ''
        go_description = ''
        pfam = ''
        kegg_pathway = ''
        cog_category = ''
        kegg_ko = ''

        # Try TCONS first in emapper, then LOC
        emap_rec = emapper.get(isoform_id) or emapper.get(ref_gene)
        if emap_rec:
            n_annot_found += 1
            raw_desc = emap_rec.get('Description', '')
            annotation = urllib.parse.unquote(raw_desc) if raw_desc else ''

            pfam = emap_rec.get('PFAMs', '')
            if pfam == '-':
                pfam = ''

            gos_raw = emap_rec.get('GOs', '')
            go_description = _resolve_top_go(gos_raw, go_map, max_terms=5)

            # KEGG pathways (map#### → readable names)
            kw_raw = emap_rec.get('KEGG_Pathway', '')
            if kw_raw and kw_raw != '-':
                kegg_pathway = _resolve_kegg_modules(kw_raw, go_map)  # reuse go_map's KEGG entries if present

            # COG category
            cog_raw = emap_rec.get('COG_category', '')
            if cog_raw and cog_raw != '-':
                cog_category = _expand_cog_letters(cog_raw)

            # KEGG KO
            ko_raw = emap_rec.get('KEGG_ko', '')
            if ko_raw and ko_raw != '-':
                kegg_ko = ko_raw
        else:
            # Fallback: try GFF description
            if ref_gene and ref_gene in gff_gene_map:
                annotation = gff_gene_map[ref_gene].get('description', '')
                if annotation:
                    annotation = urllib.parse.unquote(annotation)

        rows.append({
            'Metabolites': metabolite_name,
            'Chr': chr_name,
            'Gene_start': gene_start,
            'Gene_end': gene_end,
            'ref_Gene': ref_gene,
            'Isoform': isoform_id,
            'Pvalue': pval,
            'Annotation': annotation,
            'GO_description': go_description,
            'PFAM': pfam,
            'KEGG_Pathway': kegg_pathway,
            'COG_category': cog_category,
            'KEGG_ko': kegg_ko,
        })

    print(f"  [{metabolite_name}] {len(df):,} significant hits (p≤{pvalue_threshold})"
          f" → {len(rows)} rows  |  TCONS mapped: {n_tcons_found}"
          f"  LOC mapped: {n_loc_found}  |  Annotated: {n_annot_found}")
    return rows


def _expand_cog_letters(code):
    """Expand COG category letter(s) into human-readable names.

    ``"G"`` → ``"Carbohydrate transport and metabolism"``;
    ``"GMW"`` → ``"[G] Carbohydrate transport and metabolism; [M] Cell wall...; [W] Extracellular structures"``.
    """
    if not code or code == '-':
        return ''
    COG = {
        'A': 'RNA processing and modification',
        'B': 'Chromatin structure and dynamics',
        'C': 'Energy production and conversion',
        'D': 'Cell cycle control, cell division, chromosome partitioning',
        'E': 'Amino acid transport and metabolism',
        'F': 'Nucleotide transport and metabolism',
        'G': 'Carbohydrate transport and metabolism',
        'H': 'Coenzyme transport and metabolism',
        'I': 'Lipid transport and metabolism',
        'J': 'Translation, ribosomal structure and biogenesis',
        'K': 'Transcription',
        'L': 'Replication, recombination and repair',
        'M': 'Cell wall/membrane/envelope biogenesis',
        'N': 'Cell motility',
        'O': 'Posttranslational modification, protein turnover, chaperones',
        'P': 'Inorganic ion transport and metabolism',
        'Q': 'Secondary metabolites biosynthesis, transport and catabolism',
        'R': 'General function prediction only',
        'S': 'Function unknown',
        'T': 'Signal transduction mechanisms',
        'U': 'Intracellular trafficking, secretion, and vesicular transport',
        'V': 'Defense mechanisms',
        'W': 'Extracellular structures',
        'X': 'Mobilome: prophages, transposons',
        'Y': 'Nuclear structure',
        'Z': 'Cytoskeleton',
    }
    expanded = []
    for ch in code:
        name = COG.get(ch)
        if name:
            expanded.append(f"[{ch}] {name}")
    return '; '.join(expanded) if expanded else code


def _resolve_kegg_modules(kegg_pathway_raw, go_map=None):
    """Convert KEGG module IDs (map#####) to readable labels using go.tab.

    The eggNOG go.tab also contains KEGG entries prefixed with ``map``.
    E.g. ``M00001`` → ``Glycolysis (Embden-Meyerhof pathway), glucose => pyruvate``.
    """
    if not kegg_pathway_raw or kegg_pathway_raw == '-':
        return ''
    cleaned = []
    for term in kegg_pathway_raw.split(','):
        term = term.strip()
        if not term or term == '-':
            continue
        # Map common KEGG module IDs
        known_kegg = {
            'M00001': 'Glycolysis (Embden-Meyerhof pathway)',
            'M00002': 'Glycolysis, core module',
            'M00003': 'Gluconeogenesis',
            'M00004': 'Pentose phosphate pathway',
            'M00005': 'PRPP biosynthesis',
            'M00006': 'Pentose phosphate pathway, oxidative phase',
            'M00007': 'Pentose phosphate pathway, non-oxidative phase',
            'M00008': 'Entner-Doudoroff pathway',
            'M00009': 'Citrate cycle (TCA cycle)',
            'M00010': 'Citrate cycle, first carbon oxidation',
            'M00011': 'Citrate cycle, second carbon oxidation',
            'M00012': 'Glyoxylate cycle',
            'M00013': 'Malonate semialdehyde pathway',
            'M00014': 'Glucuronate pathway',
            'M00015': 'Proline biosynthesis',
            'M00016': 'Lysine biosynthesis',
            'M00017': 'Methionine biosynthesis',
            'M00018': 'Threonine biosynthesis',
            'M00019': 'Valine/isoleucine biosynthesis',
            'M00020': 'Serine biosynthesis',
            'M00021': 'Cysteine biosynthesis',
            'M00022': 'Shikimate pathway',
            'M00023': 'Tryptophan biosynthesis',
            'M00024': 'Phenylalanine biosynthesis',
            'M00025': 'Tyrosine biosynthesis',
            'M00026': 'Histidine biosynthesis',
            'M00027': 'GABA shunt',
            'M00028': 'Ornithine biosynthesis',
            'M00029': 'Urea cycle',
            'M00030': 'Lysine degradation',
            'M00031': 'Lysine degradation II',
            'M00032': 'Lysine degradation III',
            'M00033': 'Ectoine biosynthesis',
            'M00034': 'Methionine degradation',
            'M00035': 'Methionine salvage pathway',
            'M00036': 'Leucine degradation',
            'M00037': 'Melatonin biosynthesis',
            'M00038': 'Tryptophan degradation',
            'M00039': 'Monoterpenoid biosynthesis',
            'M00040': 'Tyrosine degradation',
            'M00041': 'β-Alanine biosynthesis',
            'M00042': 'Catecholamine biosynthesis',
            'M00043': 'Thyroid hormone biosynthesis',
            'M00044': 'Cholesterol biosynthesis',
            'M00045': 'Histidine degradation',
            'M00046': 'Pyrimidine degradation',
            'M00047': 'Creatine pathway',
            'M00048': 'Inosine monophosphate biosynthesis',
            'M00049': 'Adenine ribonucleotide biosynthesis',
            'M00050': 'Guanine ribonucleotide biosynthesis',
            'M00051': 'Uridine monophosphate biosynthesis',
            'M00052': 'Pyrimidine ribonucleotide biosynthesis',
            'M00053': 'Pyrimidine deoxyribonucleotide biosynthesis',
            'M00054': 'Glycosaminoglycan biosynthesis',
            'M00055': 'N-glycan precursor biosynthesis',
            'M00056': 'O-glycan biosynthesis',
            'M00057': 'Glycosaminoglycan degradation',
            'M00058': 'Glycosphingolipid biosynthesis',
            'M00059': 'Glycosphingolipid biosynthesis II',
            'M00060': 'Lipopolysaccharide biosynthesis',
            'M00061': 'Lipopolysaccharide biosynthesis II',
            'M00062': 'Fatty acid biosynthesis',
            'M00063': 'Fatty acid biosynthesis II',
            'M00064': 'Fatty acid degradation',
            'M00065': 'GP1-anchor biosynthesis',
            'M00066': 'Lactosylceramide biosynthesis',
            'M00067': 'Sulfoglycolipid biosynthesis',
            'M00068': 'Glycolipid biosynthesis',
            'M00069': 'Glycosphingolipid biosynthesis III',
            'M00070': 'Glycosphingolipid biosynthesis IV',
            'M00071': 'Glycosphingolipid biosynthesis V',
            'M00072': 'N-glycan biosynthesis, complex type',
            'M00073': 'N-glycan biosynthesis, high-mannose type',
            'M00074': 'N-glycan biosynthesis, hybrid type',
            'M00075': 'N-glycan biosynthesis, bisecting type',
            'M00076': 'N-glycan biosynthesis, paucimannose type',
            'M00077': 'Chondroitin sulfate biosynthesis',
            'M00078': 'Heparan sulfate biosynthesis',
            'M00079': 'Keratan sulfate biosynthesis',
            'M00080': 'Hyaluronan biosynthesis',
            'M00081': 'Sphingolipid biosynthesis',
            'M00082': 'Fatty acid biosynthesis, initiation',
            'M00083': 'Fatty acid biosynthesis, elongation',
            'M00084': 'Fatty acid biosynthesis, termination',
            'M00085': 'Fatty acid desaturation',
            'M00086': 'β-Oxidation, acyl-CoA synthesis',
            'M00087': 'β-Oxidation',
            'M00088': 'Ketone body biosynthesis',
            'M00089': 'Triacylglycerol biosynthesis',
            'M00090': 'Phosphatidylcholine biosynthesis',
            'M00091': 'Phosphatidylserine biosynthesis',
            'M00092': 'Phosphatidylethanolamine biosynthesis',
            'M00093': 'Phosphatidylinositol biosynthesis',
            'M00094': 'Ceramide biosynthesis',
            'M00095': 'C5 isoprenoid biosynthesis, mevalonate pathway',
            'M00096': 'C5 isoprenoid biosynthesis, non-mevalonate pathway',
            'M00097': 'beta-Carotene biosynthesis',
            'M00098': 'Acylglycerol degradation',
            'M00099': 'Sphingosine biosynthesis',
            'M00100': 'Bile acid biosynthesis',
            # KEGG pathway-level IDs
            'ko01100': 'Metabolic pathways',
            'ko01110': 'Biosynthesis of secondary metabolites',
            'ko01120': 'Microbial metabolism in diverse environments',
            'ko01200': 'Carbon metabolism',
            'ko01210': '2-Oxocarboxylic acid metabolism',
            'ko01212': 'Fatty acid metabolism',
            'ko01230': 'Biosynthesis of amino acids',
            'ko01232': 'Nucleotide metabolism',
            'ko01250': 'Biosynthesis of nucleotide sugars',
            'ko01240': 'Biosynthesis of cofactors',
            'ko03010': 'Ribosome',
            'ko03013': 'Nucleocytoplasmic transport',
            'ko03015': 'mRNA surveillance pathway',
            'ko03018': 'RNA degradation',
            'ko03020': 'RNA polymerase',
            'ko03022': 'Basal transcription factors',
            'ko03030': 'DNA replication',
            'ko03040': 'Spliceosome',
            'ko03050': 'Proteasome',
            'ko03060': 'Protein export',
            'ko03070': 'Bacterial secretion system',
            'ko04010': 'MAPK signaling pathway',
            'ko04011': 'MAPK signaling pathway - yeast',
            'ko04013': 'MAPK signaling pathway - fly',
            'ko04016': 'MAPK signaling pathway - plant',
            'ko04020': 'Calcium signaling pathway',
            'ko04070': 'Phosphatidylinositol signaling system',
            'ko04075': 'Plant hormone signal transduction',
            'ko04110': 'Cell cycle',
            'ko04111': 'Cell cycle - yeast',
            'ko04120': 'Ubiquitin mediated proteolysis',
            'ko04122': 'Sulfur relay system',
            'ko04130': 'SNARE interactions in vesicular transport',
            'ko04140': 'Autophagy - animal',
            'ko04141': 'Protein processing in endoplasmic reticulum',
            'ko04144': 'Endocytosis',
            'ko04145': 'Phagosome',
            'ko04146': 'Peroxisome',
            'ko04147': 'Exosome',
            'ko04151': 'PI3K-Akt signaling pathway',
            'ko04210': 'Apoptosis',
            'ko04510': 'Focal adhesion',
            'ko04626': 'Plant-pathogen interaction',
            'ko04712': 'Circadian rhythm - plant',
            'ko00500': 'Starch and sucrose metabolism',
            'ko00520': 'Amino sugar and nucleotide sugar metabolism',
            'ko00620': 'Pyruvate metabolism',
            'ko00630': 'Glyoxylate and dicarboxylate metabolism',
            'ko00900': 'Terpenoid backbone biosynthesis',
            'ko00940': 'Phenylpropanoid biosynthesis',
            'ko00941': 'Flavonoid biosynthesis',
            'ko00944': 'Flavone and flavonol biosynthesis',
            'map01100': 'Metabolic pathways (ref)',
            'map01110': 'Biosynthesis of secondary metabolites (ref)',
        }
        if term in known_kegg:
            cleaned.append(f"KEGG:{term}:{known_kegg[term]}")
        elif go_map and term in go_map:
            cleaned.append(f"KEGG:{term}:{go_map[term]}")
        else:
            cleaned.append(f"KEGG:{term}")
    return '; '.join(cleaned)


def _chr_num_to_lsat(chr_str):
    """Convert a numeric chromosome string to Lsat_Salinas_v11_chrN format.

    ``"1"`` → ``"Lsat_Salinas_v11_chr1"``;
    ``"10"``, ``"X"`` etc. are returned as-is or with Lsat_ prefix.
    """
    try:
        n = int(chr_str)
        return f'Lsat_Salinas_v11_chr{n}'
    except (ValueError, TypeError):
        return chr_str


def _resolve_top_go(gos_string, go_map, max_terms=5):
    """Resolve comma-separated GO IDs to human-readable descriptions.

    Returns the top ``max_terms`` descriptions joined by ``; ``.
    Skips root-level GO terms (molecular_function, biological_process,
    cellular_component).
    """
    if not gos_string or gos_string == '-':
        return ''

    # Root GO terms to skip
    ROOT_GOS = {
        'GO:0003674',  # molecular_function
        'GO:0008150',  # biological_process
        'GO:0005575',  # cellular_component
    }

    terms = [t.strip() for t in gos_string.split(',') if t.strip() and t.strip() != '-']
    descriptions = []
    for go in terms:
        if go in ROOT_GOS:
            continue
        desc = go_map.get(go, '')
        if desc:
            descriptions.append(desc)
        if len(descriptions) >= max_terms:
            break
    return '; '.join(descriptions)


# ---------------------------------------------------------------------------
# AI Enrichment (dual-provider: DeepSeek + OpenAI-compatible)
# ---------------------------------------------------------------------------

# Batch prompt: send all isoforms for one metabolite in a single call.
# The model sees the full picture and returns structured JSON.
_BATCH_PROMPT = """You are a plant molecular biologist analyzing GWAS results
linking metabolites to isoforms in lettuce (Lactuca sativa).

Metabolite: {metabolite}

For each significant isoform below, examine its molecular function
(PFAM domains, GO terms, and description) and hypothesize how it could
mechanistically link to {metabolite} metabolism.  If no plausible link
exists, say so explicitly.

Return ONLY a JSON array (no markdown, no preamble):
[
  {{
    "isoform": "<isoform_id>",
    "hypothesis": "<one sentence, ≤120 chars>",
    "confidence": "high|medium|low|none"
  }},
  ...
]

Isoforms:
{isoform_list}"""

_MAX_RETRIES = 3
_RETRY_BACKOFF = 2.0   # seconds, doubles each retry


def _is_openai_provider(api_base):
    """Return True if the API base URL is OpenAI (not DeepSeek or other).

    Used to toggle between OpenAI's ``reasoning`` parameter and DeepSeek's
    ``thinking`` parameter when ``--thinking`` is enabled.
    """
    return 'api.openai.com' in api_base


def _is_informative(row):
    """Return True if the row carries enough signal for AI enrichment.

    A row is *uninformative* (skipped) when ALL of these hold:
      - Annotation is empty, ``-``, or merely ``uncharacterized LOC...``
      - GO_description is empty
      - PFAM is empty
    """
    annot = (row.get('Annotation') or '').strip()
    go_d  = (row.get('GO_description') or '').strip()
    pfam  = (row.get('PFAM') or '').strip()

    # Strip useless annotation strings
    if annot.lower().startswith('uncharacterized loc') or annot in ('-', ''):
        annot = ''

    return bool(annot or go_d or pfam)


def _call_api(url, headers, payload, verbose=True, label=''):
    """POST to the chat-completions endpoint with retry + backoff.

    Returns the parsed JSON body on success, or raises on failure.
    Retries on: 429 (rate-limit), 5xx (server errors), and timeouts.
    """
    last_error = None
    for attempt in range(1, _MAX_RETRIES + 1):
        try:
            req = urllib.request.Request(
                url,
                data=json.dumps(payload).encode('utf-8'),
                headers=headers,
                method='POST',
            )
            resp = urllib.request.urlopen(req, timeout=120)
            body = json.loads(resp.read().decode('utf-8'))
            # Warn if model hit token limit (common with thinking mode)
            finish = (body.get('choices', [{}])[0]
                          .get('finish_reason', ''))
            if finish == 'length' and verbose:
                usage = body.get('usage', {})
                comp = usage.get('completion_tokens', '?')
                reasoning = (usage.get('completion_tokens_details', {})
                                  .get('reasoning_tokens', 0))
                print(f"  [{label}] ⚠  finish_reason=length "
                      f"(used {comp} completion, {reasoning} reasoning) — "
                      f"increase max_tokens to avoid truncated output")
            return body
        except urllib.error.HTTPError as e:
            status = e.code
            if status == 429 or status >= 500:
                last_error = e
                delay = _RETRY_BACKOFF ** attempt
                if verbose:
                    print(f"  [{label}] HTTP {status}, retry {attempt}/{_MAX_RETRIES}"
                          f" after {delay:.1f}s")
                time.sleep(delay)
            else:
                raise
        except (urllib.error.URLError, OSError) as e:
            last_error = e
            delay = _RETRY_BACKOFF ** attempt
            if verbose and attempt < _MAX_RETRIES:
                print(f"  [{label}] {type(e).__name__}, retry {attempt}/{_MAX_RETRIES}"
                      f" after {delay:.1f}s")
                time.sleep(delay)
    raise last_error


def _format_isoform_entry(row):
    """Format one isoform row into a compact text block for the prompt."""
    parts = [f"  - {row['Isoform']}"]
    gene = row.get('ref_Gene', '')
    if gene:
        parts.append(f"    gene: {gene}")
    pfam = (row.get('PFAM') or '').strip()
    if pfam:
        parts.append(f"    PFAM: {pfam}")
    go_desc = (row.get('GO_description') or '').strip()
    if go_desc:
        parts.append(f"    GO: {go_desc}")
    annot = (row.get('Annotation') or '').strip()
    if annot:
        parts.append(f"    desc: {annot}")
    return '\n'.join(parts)


def enrich_with_ai(rows, api_key, api_base, model, system_instruction=None,
                   batch_delay=0.1, thinking=False, reasoning_effort='medium',
                   verbose=True):
    """Enrich annotated rows with AI-generated phenotype–isoform links.

    **Batch-per-metabolite**: all isoforms for one metabolite are sent in a
    single API call.  The model returns structured JSON with a hypothesis and
    confidence level for each isoform.

    Rows added/updated keys:
      ``AI_link``       — one-sentence hypothesis (or ``-`` if skipped)
      ``AI_confidence`` — ``high``, ``medium``, ``low``, ``none``, ``-``

    Parameters
    ----------
    rows : list of dict
        All rows (both annotated and bare) from the pipeline.
    api_key : str
    api_base : str
    model : str
    batch_delay : float
        Seconds to sleep between metabolite-batch API calls.
    verbose : bool

    Returns
    -------
    list of dict
        Input rows with ``AI_link`` / ``AI_confidence`` keys added.
    """
    url = api_base.rstrip('/') + '/chat/completions'
    headers = {
        'Content-Type': 'application/json',
        'Authorization': f'Bearer {api_key}',
    }

    # Split: informative rows get enriched; bare rows get '-' default
    informative = []
    bare = []
    for r in rows:
        if _is_informative(r):
            informative.append(r)
        else:
            r['AI_link'] = '-'
            r['AI_confidence'] = '-'
            bare.append(r)

    # Group informative rows by metabolite
    from collections import defaultdict
    metabolite_groups = defaultdict(list)
    for r in informative:
        metabolite_groups[r['Metabolites']].append(r)

    n_batches = len(metabolite_groups)
    n_inf = len(informative)
    n_bare = len(bare)
    print(f"\n[AI] {n_inf} informative rows grouped into "
          f"{n_batches} metabolite batches  |  {n_bare} bare rows skipped")
    print(f"[AI]   model: {model}  @  {api_base}")

    batch_i = 0
    for metabolite, group in sorted(metabolite_groups.items()):
        batch_i += 1
        isoform_block = '\n'.join(_format_isoform_entry(r) for r in group)
        prompt = _BATCH_PROMPT.format(
            metabolite=metabolite,
            isoform_list=isoform_block,
        )

        system_msg = system_instruction or (
            'You are a concise plant molecular biologist. '
            'Return ONLY valid JSON — no markdown, no explanation.'
        )

        payload = {
            'model': model,
            'messages': [
                {'role': 'system', 'content': system_msg},
                {'role': 'user', 'content': prompt},
            ],
            'max_tokens': 1200,
            'temperature': 0.3,
        }
        if thinking:
            if _is_openai_provider(api_base):
                # OpenAI: reasoning.effort controls depth
                payload['reasoning'] = {'effort': reasoning_effort}
            else:
                # DeepSeek: thinking {type: enabled}
                payload['thinking'] = {'type': 'enabled'}
            payload['max_tokens'] = 4000  # thinking/reasoning + answer share budget
        else:
            payload['response_format'] = {'type': 'json_object'}

        label = f"{metabolite} ({batch_i}/{n_batches})"
        if verbose:
            print(f"  [{label}] {len(group)} isoforms → querying ...")

        try:
            body = _call_api(url, headers, payload, verbose=verbose, label=label)
            raw = body['choices'][0]['message']['content'].strip()

            # Parse JSON — handle markdown-wrapped or bare responses
            parsed = _parse_json_response(raw, label)
            if parsed is None:
                # Fallback: mark all as error
                for r in group:
                    r['AI_link'] = f'[parse error]'
                    r['AI_confidence'] = '-'
            else:
                # _parse_json_response returns list or dict; extract list if needed
                if isinstance(parsed, dict):
                    parsed = parsed.get('modules', parsed.get('isoforms',
                                parsed.get('results', parsed.get('data',
                                list(parsed.values())[0] if parsed else []))))
                if not isinstance(parsed, list):
                    parsed = []
                result_map = {entry.get('isoform', ''): entry for entry in parsed}
                for r in group:
                    entry = result_map.get(r['Isoform'])
                    if entry:
                        r['AI_link'] = entry.get('hypothesis', '')
                        r['AI_confidence'] = entry.get('confidence', '')
                    else:
                        r['AI_link'] = '[not in response]'
                        r['AI_confidence'] = '-'

        except Exception as e:
            if verbose:
                print(f"  [{label}] FAILED after {_MAX_RETRIES} retries: {e}")
            for r in group:
                r['AI_link'] = f'[API error: {e}]'
                r['AI_confidence'] = '-'

        if verbose:
            print(f"  [{label}] done")

        if batch_delay > 0:
            time.sleep(batch_delay)

    # Combine
    all_out = informative + bare
    n_hyp = sum(1 for r in all_out if r.get('AI_link', '-') not in ('-', '')
                and not (r.get('AI_link', '') or '').startswith('['))
    print(f"[AI]   {n_hyp}/{len(all_out)} rows with AI hypotheses")
    return all_out


def _parse_json_response(raw_text, label):
    """Robust JSON extraction from model response.

    Handles:
      - Bare JSON array/object
      - JSON inside ``` fences (markdown)
      - JSON with preamble text before the braces
      - Nested JSON objects (balanced-brace extraction, not regex)
      - Dict responses (returned as-is for per-metabolite prompts)
    """
    text = raw_text.strip()

    # Strip markdown code fences
    if text.startswith('```'):
        text = re.sub(r'^```(?:json)?\s*\n?', '', text)
        text = re.sub(r'\n?```\s*$', '', text)

    try:
        parsed = json.loads(text)
    except json.JSONDecodeError:
        # Use balanced-brace extraction: find first [{, walk until balanced
        parsed = _extract_balanced_json(text)
        if parsed is None:
            # Last resort: regex (may fail on deep nesting)
            m = re.search(r'(\[.*\]|\{.*\})', text, re.DOTALL)
            if m:
                try:
                    parsed = json.loads(m.group(1))
                except json.JSONDecodeError:
                    _print_parse_failure(text, label)
                    return None
            else:
                _print_parse_failure(text, label)
                return None

    # Handle JSON object wrapper: {"isoforms": [...]}, {"modules": [...]}, etc.
    # BUT preserve per-metabolite responses that have metadata keys
    #   (metabolite, network_summary, key_regulatory_hubs).
    if isinstance(parsed, dict):
        is_per_metab = any(k in parsed for k in
                           ('metabolite', 'network_summary', 'key_regulatory_hubs',
                            'regulatory_cascades', 'bridge_isoforms',
                            'most_pleiotropic_module'))
        if not is_per_metab:
            for key in ('isoforms', 'modules', 'results', 'data', 'hypotheses'):
                if key in parsed:
                    parsed = parsed[key]
                    break
    if not isinstance(parsed, list):
        # Dict response (per-metabolite or network prompt) — return as-is
        if isinstance(parsed, dict):
            return parsed
        print(f"  [{label}] WARNING: expected JSON array or dict, got {type(parsed).__name__}")
        return None

    return parsed


def _extract_balanced_json(text):
    """Find the first `{` or `[` and extract until braces are balanced.

    Returns the parsed JSON object, or None on failure.
    """
    for start_char in ('{', '['):
        pos = text.find(start_char)
        if pos == -1:
            continue
        close_char = '}' if start_char == '{' else ']'
        depth = 0
        in_string = False
        escape = False
        for i in range(pos, len(text)):
            ch = text[i]
            if escape:
                escape = False
                continue
            if ch == '\\' and in_string:
                escape = True
                continue
            if ch == '"':
                in_string = not in_string
                continue
            if in_string:
                continue
            if ch == start_char:
                depth += 1
            elif ch == close_char:
                depth -= 1
                if depth == 0:
                    try:
                        return json.loads(text[pos:i + 1])
                    except json.JSONDecodeError:
                        return None
        # Unbalanced — try next start char
    return None


def _print_parse_failure(text, label):
    """Print diagnostic info when JSON parsing fails."""
    print(f"  [{label}] WARNING: could not parse JSON from response")
    preview = text[:200].replace('\n', '\\n')
    print(f"  [{label}]   raw preview: {preview}")


# ---------------------------------------------------------------------------
# Co-function Network + AI Regulatory Analysis
# ---------------------------------------------------------------------------

_NETWORK_ANALYSIS_PROMPT = """You are a plant systems biologist analyzing a
co-function network built from GWAS-significant isoforms across multiple
metabolite phenotypes in lettuce (Lactuca sativa).

Nodes are isoforms (transcripts/genes) that were significantly associated with
one or more metabolites.  Edges connect isoforms that share GO terms, PFAM
domains, or other functional annotations — indicating they participate in
related biological processes.

Below are the detected functional modules (connected components in the
co-function network).  For each module:

1. Identify the dominant biological function from the shared GO/PFAM/COG terms.
2. Based on the associated metabolites, hypothesize the regulatory role this
   module plays — e.g. which metabolic pathway it controls, whether it is a
   hub for multiple metabolites, or a specific regulator for one.
3. Name the module with a concise functional label (e.g. "Amino acid
   biosynthesis hub", "Sugar transport & signaling", "Lipid modification
   cascade").
4. Assign a confidence level: high (clear mechanistic link), medium (plausible
   but indirect), low (speculative).

Network-level analysis:
- Which modules likely form a regulatory cascade (one module upstream of
  another)?
- Are there "bridge" isoforms that connect otherwise-separate modules? List
  any bridge isoforms by ID.
- Which module is the most pleiotropic (associated with the most metabolites)?

Return ONLY a JSON object (no markdown, no preamble):
{
  "network_summary": "<one paragraph overview of the regulatory landscape>",
  "modules": [
    {
      "module_id": "M1",
      "name": "<functional label>",
      "size": <number of isoforms>,
      "metabolites": ["<metabolite1>", ...],
      "dominant_function": "<GO/PFAM-based function>",
      "regulatory_hypothesis": "<1-2 sentence mechanistic hypothesis>",
      "confidence": "high|medium|low",
      "hub_isoforms": ["<isoform_id>", ...],
      "potential_targets": ["<metabolic pathway or process>", ...]
    }
  ],
  "regulatory_cascades": [
    {"upstream_module": "M1", "downstream_module": "M2",
     "mechanism": "<how they are linked>"}
  ],
  "bridge_isoforms": [
    {"isoform": "<id>", "connects_modules": ["M1", "M2"],
     "function": "<brief description>"}
  ],
  "most_pleiotropic_module": "M1",
  "key_regulatory_hubs": ["<isoform_id>", ...]
}

IMPORTANT: Replace ALL placeholder values (like "M1", "<isoform_id>", etc.)
with actual module IDs and isoform IDs from the network data.  For
most_pleiotropic_module, pick the module whose metabolite list is the
longest — count the metabolites yourself and choose the real module ID.

Network data:
{module_descriptions}"""


# Generic GO terms filtered from edge-building AND module characterization
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

# COG categories too broad for co-function edge-building — they are biologically
# meaningful as annotations but spanning thousands of genes, so sharing the same
# COG letter does NOT imply a specific co-function relationship.
#
# These are the COG expanded strings as produced by _expand_cog_letters().
_GENERIC_COG_EDGE_TERMS = {
    '[K] Transcription',                         # covers all TFs — too broad
    #'[T] Signal transduction mechanisms',        # covers all kinases/receptors
    '[R] General function prediction only',      # literally generic
    '[S] Function unknown',                      # unknown function
}

# KEGG pathway IDs too broad for co-function edge-building — reference maps and
# top-level overview pathways that are assigned to nearly every metabolic gene.
_BROAD_KEGG_IDS = {
    'map01100', 'ko01100',     # Metabolic pathways (global overview)
    'map01110', 'ko01110',     # Biosynthesis of secondary metabolites (broad)
}


def _is_kegg_term_broad(term):
    """Return True if *term* is a broad KEGG reference-map entry.

    Checks the embedded pathway ID (e.g. ``map01100``) against *_BROAD_KEGG_IDS*.
    The term format is ``KEGG:map01100:Metabolic pathways (ref)``.
    """
    # Extract the path:ID portion — everything between "KEGG:" and the next ":"
    if not term.startswith('KEGG:'):
        return False
    rest = term[5:]                     # strip "KEGG:"
    colon = rest.find(':')
    if colon == -1:
        # No description suffix, e.g. "KEGG:map01100"
        return rest in _BROAD_KEGG_IDS
    path_id = rest[:colon]              # e.g. "map01100"
    return path_id in _BROAD_KEGG_IDS


def build_cofunction_network(all_rows, top_n=300, min_shared=1):
    """Build a co-function network from annotated GWAS hits.

    Isoforms are connected when they share GO terms or PFAM domains — a
    stronger signal than sharing a broad COG category.  Edges are weighted
    by the number of shared terms.

    Parameters
    ----------
    all_rows : list of dict
        Annotated rows from the pipeline.
    top_n : int
        Take the top-N most significant isoforms (by p-value).
    min_shared : int
        Minimum shared terms to create an edge.

    Returns
    -------
    nodes : dict
        isoform_id → {Isoform, ref_Gene, Pvalue, Annotation, PFAM, GO_description,
                       Metabolites, Chr, Gene_start, Gene_end}
    edges : list of dict
        {source, target, weight, shared_terms, term_type}
    modules : list of list
        Each inner list is a connected component (list of isoform IDs).
    metabolite_map : dict
        isoform_id → set of associated metabolite names.
    """
    print(f"\n[NETWORK] Building co-function network from top {top_n} isoforms …")

    # Sort by p-value, take top-N, deduplicate by isoform
    seen = set()
    top_rows = []
    for r in sorted(all_rows, key=lambda x: x.get('Pvalue', 1.0)):
        iso = r['Isoform']
        if iso not in seen:
            seen.add(iso)
            top_rows.append(r)
        if len(top_rows) >= top_n:
            break

    # Index by isoform ID
    nodes = {}
    for r in top_rows:
        nodes[r['Isoform']] = {
            'Isoform': r['Isoform'],
            'ref_Gene': r.get('ref_Gene', ''),
            'Pvalue': r['Pvalue'],
            'Annotation': r.get('Annotation', ''),
            'PFAM': r.get('PFAM', ''),
            'GO_description': r.get('GO_description', ''),
            'KEGG_Pathway': r.get('KEGG_Pathway', ''),
            'COG_category': r.get('COG_category', ''),
            'KEGG_ko': r.get('KEGG_ko', ''),
            'Metabolites': r.get('Metabolites', ''),
            'Chr': r.get('Chr', ''),
            'Gene_start': r.get('Gene_start', ''),
            'Gene_end': r.get('Gene_end', ''),
        }

    # Build term → isoforms index for ALL annotation types
    go_index = {}      # term → set of isoform IDs
    pfam_index = {}    # term → set of isoform IDs
    kegg_index = {}    # term → set of isoform IDs
    cog_index = {}     # term → set of isoform IDs
    desc_index = {}    # term → set of isoform IDs (non-generic descriptions only)

    # Generic description patterns to skip (too vague for co-function signal)
    _VAGUE_DESC_RE = re.compile(
        r'^(uncharacteri[sz]ed|hypothetical|predicted|probable|putative|'
        r'expressed|unknown|unnamed|orf|cds|partial|'
        r'encoded by|belongs to the|contain|domain|family|'
        r'protein of unknown|duf\d+)',
        re.IGNORECASE)

    for iso_id, nd in nodes.items():
        # GO terms (from GO_description — semicolon-separated)
        go_text = nd.get('GO_description', '')
        if go_text:
            for term in go_text.split(';'):
                term = term.strip().lower()
                if term and term != '-':
                    go_index.setdefault(term, set()).add(iso_id)

        # PFAM domains (comma-separated)
        pfam_text = nd.get('PFAM', '')
        if pfam_text:
            for term in pfam_text.split(','):
                term = term.strip()
                if term and term != '-':
                    pfam_index.setdefault(term, set()).add(iso_id)

        # KEGG pathway modules (semicolon-separated: "KEGG:M00001:Glycolysis; ...")
        kegg_text = nd.get('KEGG_Pathway', '')
        if kegg_text:
            for term in kegg_text.split(';'):
                term = term.strip()
                if term and term != '-':
                    kegg_index.setdefault(term, set()).add(iso_id)

        # COG categories (semicolon-separated: "[G] Carbohydrate...; [E] Amino acid...")
        cog_text = nd.get('COG_category', '')
        if cog_text:
            for term in cog_text.split(';'):
                term = term.strip()
                if term and term != '-':
                    cog_index.setdefault(term, set()).add(iso_id)

        # eggNOG Description — use as term ONLY if specific (not vague)
        desc_text = nd.get('Annotation', '')
        if desc_text and not _VAGUE_DESC_RE.search(desc_text):
            # Truncate long descriptions for indexing
            desc_key = desc_text[:80].strip().lower()
            if desc_key:
                desc_index.setdefault(desc_key, set()).add(iso_id)

    n_go_removed = 0
    for term in list(go_index.keys()):
        if term in _GENERIC_GO_TERMS:
            del go_index[term]
            n_go_removed += 1

    # Filter broad COG categories from edge-building (same rationale as GO)
    n_cog_removed = 0
    for term in list(cog_index.keys()):
        if term in _GENERIC_COG_EDGE_TERMS:
            del cog_index[term]
            n_cog_removed += 1

    # Filter broad KEGG reference maps from edge-building
    n_kegg_removed = 0
    for term in list(kegg_index.keys()):
        if _is_kegg_term_broad(term):
            del kegg_index[term]
            n_kegg_removed += 1

    # Build edges: isoform pairs sharing terms
    edge_weights = {}  # (iso_a, iso_b) → {weight, shared_terms, term_type}

    for term_index, term_type in [(go_index, 'GO'), (pfam_index, 'PFAM'),
                                   (kegg_index, 'KEGG'), (cog_index, 'COG'),
                                   (desc_index, 'DESC')]:
        for term, iso_set in term_index.items():
            iso_list = sorted(iso_set)
            for i in range(len(iso_list)):
                for j in range(i + 1, len(iso_list)):
                    a, b = iso_list[i], iso_list[j]
                    key = (a, b) if a < b else (b, a)
                    if key not in edge_weights:
                        edge_weights[key] = {'weight': 0, 'shared_terms': [],
                                              'term_type': set()}
                    edge_weights[key]['weight'] += 1
                    edge_weights[key]['shared_terms'].append(f"{term_type}:{term}")
                    edge_weights[key]['term_type'].add(term_type)

    edges = []
    for (a, b), ew in edge_weights.items():
        if ew['weight'] >= min_shared:
            edges.append({
                'source': a,
                'target': b,
                'weight': ew['weight'],
                'shared_terms': '; '.join(ew['shared_terms']),
                'term_type': ','.join(sorted(ew['term_type'])),
            })

    # Build adjacency for connected components
    adj = {}
    for iso_id in nodes:
        adj[iso_id] = []
    for e in edges:
        adj[e['source']].append(e['target'])
        adj[e['target']].append(e['source'])

    # Connected components (modules)
    visited = set()
    modules = []
    for iso_id in nodes:
        if iso_id not in visited:
            comp = []
            stack = [iso_id]
            while stack:
                v = stack.pop()
                if v not in visited:
                    visited.add(v)
                    comp.append(v)
                    for nb in adj.get(v, []):
                        if nb not in visited:
                            stack.append(nb)
            modules.append(comp)

    # Sort modules by size descending
    modules.sort(key=len, reverse=True)

    # Build metabolite map
    metabolite_map = {}
    for iso_id, nd in nodes.items():
        mets = set(m.strip() for m in nd['Metabolites'].split(','))
        metabolite_map[iso_id] = mets

    n_edges_total = len(edges)
    n_solo = sum(1 for m in modules if len(m) == 1)
    n_mod = sum(1 for m in modules if len(m) > 1)
    avg_size = sum(len(m) for m in modules) / max(len(modules), 1)

    print(f"[NETWORK]   {len(nodes)} nodes  |  {n_edges_total} edges"
          f"  |  {n_mod} modules (≥2 nodes)  |  {n_solo} singletons"
          f"  |  avg module size: {avg_size:.1f}")

    return nodes, edges, modules, metabolite_map


def _describe_module(module_isoforms, nodes, metabolite_map, module_id):
    """Build a text description of one network module for the AI prompt."""
    lines = [f"Module {module_id} ({len(module_isoforms)} isoforms):"]

    # Collect metabolites
    all_mets = set()
    for iso in module_isoforms:
        all_mets |= metabolite_map.get(iso, set())
    lines.append(f"  Associated metabolites: {', '.join(sorted(all_mets)[:15])}")

    # Collect top functional terms
    go_terms = []
    pfam_terms = []
    for iso in module_isoforms:
        nd = nodes[iso]
        go = nd.get('GO_description', '')
        if go:
            for t in go.split(';'):
                t = t.strip()
                if t and t != '-':
                    go_terms.append(t)
        pf = nd.get('PFAM', '')
        if pf:
            for t in pf.split(','):
                t = t.strip()
                if t and t != '-':
                    pfam_terms.append(t)

    from collections import Counter
    top_go = [t for t, _ in Counter(go_terms).most_common(5)]
    top_pfam = [t for t, _ in Counter(pfam_terms).most_common(5)]

    lines.append(f"  Top GO terms: {', '.join(top_go) if top_go else '(none)'}")
    lines.append(f"  Top PFAM domains: {', '.join(top_pfam) if top_pfam else '(none)'}")

    # List isoforms with annotations
    lines.append("  Isoforms:")
    for iso in module_isoforms:
        nd = nodes[iso]
        annot = nd.get('Annotation', '')[:60]
        desc = f"{iso} (gene: {nd.get('ref_Gene','?')}, p={nd['Pvalue']:.2e})"
        if annot:
            desc += f" — {annot}"
        lines.append(f"    {desc}")

    return '\n'.join(lines)


def analyze_network_with_ai(nodes, edges, modules, metabolite_map,
                             api_key, api_base, model,
                             min_module_size=3, batch_delay=0.1,
                             thinking=False, reasoning_effort='medium',
                             verbose=True):
    """Send network module descriptions to AI for regulatory interpretation.

    Parameters
    ----------
    nodes, edges, modules, metabolite_map :
        Outputs of ``build_cofunction_network``.
    api_key, api_base, model : str
        API configuration.
    min_module_size : int
        Modules smaller than this are only listed, not analyzed in depth.
    batch_delay : float
        Seconds to sleep before the API call.
    verbose : bool

    Returns
    -------
    dict
        Parsed AI response with keys: network_summary, modules, regulatory_cascades,
        bridge_isoforms, most_pleiotropic_module, key_regulatory_hubs.
        Returns ``None`` on failure.
    """
    big_modules = [(i, m) for i, m in enumerate(modules, 1)
                   if len(m) >= min_module_size]
    small_modules = [(i, m) for i, m in enumerate(modules, 1)
                     if len(m) < min_module_size]

    if not big_modules:
        print("[NETWORK-AI] No modules meet minimum size — skipping AI analysis.")
        return None

    # Build prompt descriptions for big modules
    descriptions = []
    for mod_id, mod_iso in big_modules:
        descriptions.append(_describe_module(mod_iso, nodes, metabolite_map,
                                              f"M{mod_id}"))
    # Add a summary line for small modules
    if small_modules:
        small_summary = "Small modules (not analyzed in depth):\n"
        for mod_id, mod_iso in small_modules:
            all_mets = set()
            for iso in mod_iso:
                all_mets |= metabolite_map.get(iso, set())
            small_summary += (f"  M{mod_id}: {len(mod_iso)} isoforms, "
                              f"metabolites: {', '.join(sorted(all_mets)[:8])}\n")
        descriptions.append(small_summary)

    module_block = '\n\n'.join(descriptions)
    prompt = _NETWORK_ANALYSIS_PROMPT.format(module_descriptions=module_block)

    url = api_base.rstrip('/') + '/chat/completions'
    headers = {
        'Content-Type': 'application/json',
        'Authorization': f'Bearer {api_key}',
    }

    payload = {
        'model': model,
        'messages': [
            {'role': 'system',
             'content': 'You are a plant systems biologist. '
                        'Return ONLY valid JSON — no markdown, no explanation.'},
            {'role': 'user', 'content': prompt},
        ],
        'max_tokens': 4000,
        'temperature': 0.3,
        'response_format': {'type': 'json_object'},
    }
    if thinking:
        if _is_openai_provider(api_base):
            payload['reasoning'] = {'effort': reasoning_effort}
        else:
            payload['thinking'] = {'type': 'enabled'}
        payload['max_tokens'] = 16000  # thinking/reasoning + answer share budget
        payload.pop('response_format', None)

    n_mod_big = len(big_modules)
    print(f"\n[NETWORK-AI] Sending {n_mod_big} modules ({sum(len(m) for _, m in big_modules)}"
          f" isoforms) to {model} for regulatory analysis …")

    try:
        body = _call_api(url, headers, payload, verbose=verbose,
                         label=f"network ({n_mod_big} modules)")
        raw = body['choices'][0]['message']['content'].strip()
        parsed = json.loads(raw)
        print(f"[NETWORK-AI]   analysis complete")
        return parsed
    except json.JSONDecodeError:
        # Try _parse_json_response for markdown-wrapped JSON
        parsed = _parse_json_response(raw, f"network ({n_mod_big} modules)")
        if parsed:
            if isinstance(parsed, list):
                parsed = {'modules': parsed}
            print(f"[NETWORK-AI]   analysis complete (fallback parse)")
            return parsed
        print(f"[NETWORK-AI]   ERROR: could not parse AI response as JSON")
        return None
    except Exception as e:
        print(f"[NETWORK-AI]   ERROR: {e}")
        return None


def _sanitize_ai_response(ai_analysis, modules, metabolite_map):
    """Fix placeholder values the model may have copied verbatim from the prompt.

    Detects and replaces:
      - ``most_pleiotropic_module`` values of ``"M?"``, ``"?"``, or missing
      - Any other obviously-stale template values.
    """
    _PLACEHOLDER_VALS = {'M?', '?', '', None}

    # --- most_pleiotropic_module ---
    pleio = ai_analysis.get('most_pleiotropic_module')
    if pleio in _PLACEHOLDER_VALS:
        # Compute from actual module data: find the module with the most
        # metabolites (by unique metabolite count).
        best_id = '?'
        best_n = 0
        for i, mod_iso in enumerate(modules, 1):
            mets = set()
            for iso in mod_iso:
                mets |= metabolite_map.get(iso, set())
            if len(mets) > best_n:
                best_n = len(mets)
                best_id = f'M{i}'
        if best_n > 0:
            ai_analysis['most_pleiotropic_module'] = best_id
        # else leave as-is (edge case: no modules at all)

    return ai_analysis


def export_network(nodes, edges, modules, metabolite_map,
                   ai_analysis, out_prefix):
    """Write all co-function network output files.

    Parameters
    ----------
    nodes, edges, modules, metabolite_map :
        Outputs of ``build_cofunction_network``.
    ai_analysis : dict or None
        Output of ``analyze_network_with_ai``.
    out_prefix : str
        Output file prefix.
    """
    print(f"\n[NETWORK] Exporting co-function network files …")

    # --- Nodes table ---
    nodes_path = f"{out_prefix}_cofunction_nodes.tsv"
    nodes_rows = []
    for iso_id, nd in nodes.items():
        mets = metabolite_map.get(iso_id, set())
        nodes_rows.append({
            'Isoform': iso_id,
            'ref_Gene': nd.get('ref_Gene', ''),
            'Chr': nd.get('Chr', ''),
            'Gene_start': nd.get('Gene_start', ''),
            'Gene_end': nd.get('Gene_end', ''),
            'Pvalue': nd.get('Pvalue', ''),
            'Annotation': (nd.get('Annotation', '') or '')[:200],
            'PFAM': nd.get('PFAM', ''),
            'GO_description': nd.get('GO_description', ''),
            'KEGG_Pathway': nd.get('KEGG_Pathway', ''),
            'COG_category': nd.get('COG_category', ''),
            'KEGG_ko': nd.get('KEGG_ko', ''),
            'Metabolites': ', '.join(sorted(mets)),
        })
    pd.DataFrame(nodes_rows).to_csv(nodes_path, sep='\t', index=False)
    print(f"  → {nodes_path}  ({len(nodes_rows)} nodes)")

    # --- Edges table ---
    edges_path = f"{out_prefix}_cofunction_edges.tsv"
    pd.DataFrame(edges).to_csv(edges_path, sep='\t', index=False)
    print(f"  → {edges_path}  ({len(edges)} edges)")

    # --- Module membership table (with characterizing terms) ---
    mod_path = f"{out_prefix}_cofunction_modules.tsv"
    mod_rows = []
    for i, mod_iso in enumerate(modules, 1):
        all_mets = set()
        # Collect all terms used by this module's isoforms
        term_counter = {}
        for iso in mod_iso:
            all_mets |= metabolite_map.get(iso, set())
            nd = nodes.get(iso, {})
            for src, sep in [('PFAM', ','), ('GO_description', ';'),
                              ('KEGG_Pathway', ';'), ('COG_category', ';')]:
                text = nd.get(src, '')
                if text:
                    for t in text.split(sep):
                        t = t.strip()
                        if not t or t in ('', '-'):
                            continue
                        # Skip generic GO terms
                        if src == 'GO_description' and t.lower() in _GENERIC_GO_TERMS:
                            continue
                        # Skip Function unknown COG
                        if src == 'COG_category' and 'Function unknown' in t:
                            continue
                        # Skip broad KEGG reference maps (same filter as edge-building)
                        if src == 'KEGG_Pathway' and _is_kegg_term_broad(t):
                            continue
                        tag = f"{src}:{t}"
                        term_counter[tag] = term_counter.get(tag, 0) + 1
            # Also include brief Annotation as a characterizing descriptor
            annot = nd.get('Annotation', '')
            if annot and len(annot) > 3 and len(annot) < 80 \
               and not annot.lower().startswith('uncharacterized') \
               and not annot.lower().startswith('encoded by') \
               and not annot.lower().startswith('belongs to'):
                tag = f"DESC:{annot}"
                term_counter[tag] = term_counter.get(tag, 0) + 1
        # Top terms (most isoforms share them)
        top_terms = sorted(term_counter.items(), key=lambda x: -x[1])[:8]
        top_terms_str = '; '.join(f"{t} ({c})" for t, c in top_terms)
        mod_rows.append({
            'module_id': f'M{i}',
            'size': len(mod_iso),
            'isoforms': ', '.join(mod_iso),
            'metabolites': ', '.join(sorted(all_mets)),
            'top_terms': top_terms_str,
        })
    pd.DataFrame(mod_rows).to_csv(mod_path, sep='\t', index=False)
    print(f"  → {mod_path}  ({len(mod_rows)} modules)")

    # --- AI regulatory analysis report ---
    if ai_analysis:
        # Sanitize placeholder values the model may have copied verbatim
        ai_analysis = _sanitize_ai_response(ai_analysis, modules, metabolite_map)
        report_path = f"{out_prefix}_network_ai_report.txt"
        with open(report_path, 'w') as f:
            f.write("=" * 78 + "\n")
            f.write("  Co-Function Network -- AI Regulatory Analysis\n")
            f.write("=" * 78 + "\n\n")

            # Network summary
            ns = ai_analysis.get('network_summary', '')
            if ns:
                f.write("NETWORK OVERVIEW\n")
                f.write("-" * 78 + "\n")
                f.write(ns + "\n\n")

            # Modules
            ai_mods = ai_analysis.get('modules', [])
            if ai_mods:
                f.write("FUNCTIONAL MODULES\n")
                f.write("-" * 78 + "\n\n")
                for am in ai_mods:
                    f.write(f"  {am.get('module_id', '?')}: "
                            f"{am.get('name', 'Unnamed')}\n")
                    f.write(f"    Size: {am.get('size', '?')} isoforms  |  "
                            f"Confidence: {am.get('confidence', '?')}\n")
                    mets = am.get('metabolites', [])
                    if mets:
                        f.write(f"    Metabolites: {', '.join(mets)}\n")
                    f.write(f"    Function: {am.get('dominant_function', '?')}\n")
                    f.write(f"    Hypothesis: {am.get('regulatory_hypothesis', '?')}\n")
                    hubs = am.get('hub_isoforms', [])
                    if hubs:
                        f.write(f"    Hub isoforms: {', '.join(hubs)}\n")
                    targets = am.get('potential_targets', [])
                    if targets:
                        f.write(f"    Potential targets: {', '.join(targets)}\n")
                    f.write("\n")

            # Regulatory cascades
            cascades = ai_analysis.get('regulatory_cascades', [])
            if cascades:
                f.write("REGULATORY CASCADES\n")
                f.write("-" * 78 + "\n\n")
                for c in cascades:
                    f.write(f"  {c.get('upstream_module','?')} → "
                            f"{c.get('downstream_module','?')}\n")
                    f.write(f"    Mechanism: {c.get('mechanism', '?')}\n\n")

            # Bridge isoforms
            bridges = ai_analysis.get('bridge_isoforms', [])
            if bridges:
                f.write("BRIDGE ISOFORMS (connect multiple modules)\n")
                f.write("-" * 78 + "\n\n")
                for b in bridges:
                    f.write(f"  {b.get('isoform', '?')}\n")
                    f.write(f"    Connects: {', '.join(b.get('connects_modules', []))}\n")
                    f.write(f"    Function: {b.get('function', '?')}\n\n")

            # Key findings
            f.write("KEY FINDINGS\n")
            f.write("-" * 78 + "\n\n")
            pleio = ai_analysis.get('most_pleiotropic_module', '?')
            f.write(f"  Most pleiotropic module: {pleio}\n")
            hubs = ai_analysis.get('key_regulatory_hubs', [])
            if hubs:
                f.write(f"  Key regulatory hubs: {', '.join(hubs)}\n")

        print(f"  → {report_path}")

        # JSON export for downstream tools
        json_path = f"{out_prefix}_network_ai_analysis.json"
        with open(json_path, 'w') as f:
            json.dump(ai_analysis, f, indent=2, default=str)
        print(f"  → {json_path}")


def run_network_pipeline(qqman_files, metabolite_names, args,
                         tcons_map, gff_gene_map, emapper, go_map):
    """Build one co-function network **per metabolite**.

    Each metabolite's top ``--network-top-n`` isoforms form their own
    network — isoforms that co-associate with the SAME metabolite and
    share GO/PFAM terms are linked.  This preserves pathway specificity:
    fucose isoforms connect to fucose-related functions, citrate isoforms
    to TCA-cycle functions, etc.

    A cross-metabolite summary identifies pleiotropic isoforms that appear
    in multiple metabolites' networks.
    """
    per_metab_top = args.network_per_metab_top

    # 1. Collect isoforms per metabolite (no pooling)
    metab_isoforms = {}   # metabolite → list of annotated row dicts
    for qqman_path, met_name in zip(qqman_files, metabolite_names):
        rows = _read_one_metabolite_top(qqman_path, met_name,
                                        per_metab_top, tcons_map,
                                        gff_gene_map, emapper, go_map)
        if rows:
            metab_isoforms[met_name] = rows

    n_metab_with_hits = len(metab_isoforms)
    total_isoforms = sum(len(v) for v in metab_isoforms.values())
    print(f"\n[NETWORK] {n_metab_with_hits} metabolites with isoforms "
          f"  ({total_isoforms:,} total isoform entries)")

    # 2. Build per-metabolite networks
    all_networks = {}  # metabolite → {nodes, edges, modules, metabolite_map}
    all_ai = {}        # metabolite → AI analysis dict or None
    isoform_presence = {}  # isoform_id → list of metabolites

    for met_name, rows in sorted(metab_isoforms.items()):
        # Track isoform→metabolite presence for cross-metabolite summary
        for r in rows:
            iso = r['Isoform']
            isoform_presence.setdefault(iso, []).append(met_name)

        if len(rows) < 2:
            continue  # need ≥2 isoforms to form any edge

        nodes, edges, modules, met_map = build_cofunction_network(
            rows, top_n=args.network_top_n,
            min_shared=args.network_min_shared,
        )
        all_networks[met_name] = {
            'nodes': nodes, 'edges': edges, 'modules': modules,
            'metabolite_map': met_map,
        }

    n_networks = len(all_networks)
    print(f"[NETWORK]   {n_networks} per-metabolite networks built")

    # 3. AI analysis per metabolite (only those with big enough modules)
    if args.api_key:
        ai_batch = {}
        for met_name, net in all_networks.items():
            big = [m for m in net['modules'] if len(m) >= args.network_min_module]
            if big:
                ai_batch[met_name] = net
        if ai_batch:
            print(f"\n[NETWORK-AI] Analyzing {len(ai_batch)} metabolite networks …")
            all_ai = _analyze_per_metabolite(
                ai_batch, args.api_key, args.api_base, args.model,
                min_module_size=args.network_min_module,
                batch_delay=args.ai_rate_limit,
                thinking=args.thinking,
                reasoning_effort=args.reasoning_effort,
            )

    # 4. Export per-metabolite network files
    exported_count = 0
    for met_name, net in sorted(all_networks.items()):
        ai = all_ai.get(met_name)
        modules = net['modules']
        # Only export if there's at least one meaningful module
        has_module = any(len(m) >= 2 for m in modules)
        if not has_module and not ai:
            continue
        safe_name = re.sub(r'[^a-zA-Z0-9_-]', '_', met_name)
        met_prefix = f"{args.out_prefix}_{safe_name}"
        export_network(net['nodes'], net['edges'], net['modules'],
                       net['metabolite_map'], ai, met_prefix)
        exported_count += 1

    print(f"[NETWORK]   exported {exported_count} per-metabolite networks")

    # 5. Cross-metabolite summary
    _export_cross_metabolite_summary(isoform_presence, all_networks, args.out_prefix)

    # 6. R plots — top metabolites by edge count
    if args.plot_network and all_networks:
        _plot_top_metabolite_networks(all_networks, all_ai, args)

    return all_networks, all_ai


def _read_one_metabolite_top(qqman_path, met_name, top_n,
                             tcons_map, gff_gene_map, emapper, go_map):
    """Read one metabolite QQMan, take top-N isoforms, annotate.

    Returns list of dicts in same format as ``process_qqman`` output.
    """
    df = pd.read_csv(
        qqman_path, sep='\t', header=None,
        names=['original_id', 'chr', 'pos', 'pvalue'],
        dtype={'original_id': str, 'chr': str, 'pos': int, 'pvalue': float}
    )
    df['gene_id'] = df['original_id'].apply(extract_gene_id)
    df = df.loc[df.groupby('gene_id')['pvalue'].idxmin()].copy()
    df = df.sort_values('pvalue', ascending=True).head(top_n)

    rows = []
    for _, row in df.iterrows():
        gid = row['gene_id']
        pval = row['pvalue']

        if gid.startswith('TCONS'):
            tinfo = tcons_map.get(gid, {})
            chr_name = tinfo.get('chr', _chr_num_to_lsat(row['chr']))
            ref_gene = tinfo.get('gene_name', '')
            if ref_gene and ref_gene in gff_gene_map:
                gs = gff_gene_map[ref_gene]['start']
                ge = gff_gene_map[ref_gene]['end']
            elif tinfo:
                gs, ge = tinfo.get('start', row['pos']), tinfo.get('end', row['pos'])
            else:
                gs = ge = row['pos']
        elif gid.startswith('LOC'):
            chr_name = _chr_num_to_lsat(row['chr'])
            ref_gene = gid
            if gid in gff_gene_map:
                gs = gff_gene_map[gid]['start']
                ge = gff_gene_map[gid]['end']
                chr_name = gff_gene_map[gid].get('chr', chr_name)
            else:
                gs = ge = row['pos']
        else:
            continue

        emap_rec = emapper.get(gid) or emapper.get(ref_gene, {})
        annotation = urllib.parse.unquote(emap_rec.get('Description', ''))
        pfam = emap_rec.get('PFAMs', '')
        if pfam == '-':
            pfam = ''
        go_desc = _resolve_top_go(emap_rec.get('GOs', ''), go_map)
        kegg_pathway = _resolve_kegg_modules(emap_rec.get('KEGG_Pathway', ''))
        cog_category = _expand_cog_letters(emap_rec.get('COG_category', ''))
        kegg_ko = emap_rec.get('KEGG_ko', '')
        if kegg_ko == '-':
            kegg_ko = ''

        rows.append({
            'Metabolites': met_name,
            'Chr': chr_name,
            'Gene_start': gs, 'Gene_end': ge,
            'ref_Gene': ref_gene, 'Isoform': gid,
            'Pvalue': pval,
            'Annotation': annotation,
            'GO_description': go_desc,
            'PFAM': pfam,
            'KEGG_Pathway': kegg_pathway,
            'COG_category': cog_category,
            'KEGG_ko': kegg_ko,
        })
    return rows


# Per-metabolite AI batch prompt
_PER_METAB_AI_PROMPT = """You are a plant systems biologist. For the metabolite
below, analyze the co-function network of its top GWAS-significant isoforms.

Metabolite: {metabolite}

Each isoform was significantly associated with this metabolite in GWAS.
Isoforms are connected in the network when they share specific GO terms or
PFAM domains — indicating they participate in related molecular functions
that may jointly regulate this metabolite's biosynthesis, transport, or
degradation.

Functional modules detected:
{module_descriptions}

For EACH module:
1. Name it with a concise functional label.
2. Hypothesize how this module regulates the metabolite (biosynthesis?
   transport? signaling? degradation?).
3. Which isoforms are likely hubs?
4. Assign confidence: high | medium | low.

Return ONLY a JSON object:
{{
  "metabolite": "{metabolite}",
  "network_summary": "<one paragraph: how these modules may jointly regulate
                       {metabolite}>",
  "modules": [
    {{
      "module_id": "M1",
      "name": "<functional label>",
      "size": <N>,
      "dominant_function": "<GO/PFAM-based function>",
      "regulatory_hypothesis": "<1-2 sentence mechanism>",
      "confidence": "high|medium|low",
      "hub_isoforms": ["<id>", ...]
    }}
  ],
  "key_regulatory_hubs": ["<isoform_id>", ...]
}}"""


def _analyze_per_metabolite(networks, api_key, api_base, model,
                            min_module_size=3, batch_delay=0.1,
                            thinking=False, reasoning_effort='medium'):
    """AI analysis: one API call per metabolite with meaningful modules."""
    url = api_base.rstrip('/') + '/chat/completions'
    headers = {
        'Content-Type': 'application/json',
        'Authorization': f'Bearer {api_key}',
    }
    results = {}
    total = len(networks)
    for i, (met_name, net) in enumerate(sorted(networks.items()), 1):
        big_mods = [(j, m) for j, m in enumerate(net['modules'], 1)
                     if len(m) >= min_module_size]
        if not big_mods:
            continue

        # Build module descriptions
        descriptions = []
        for mod_j, mod_iso in big_mods:
            descriptions.append(
                _describe_module(mod_iso, net['nodes'],
                                 net['metabolite_map'], f"M{mod_j}")
            )
        mod_block = '\n\n'.join(descriptions)
        prompt = _PER_METAB_AI_PROMPT.format(
            metabolite=met_name, module_descriptions=mod_block)

        payload = {
            'model': model,
            'messages': [
                {'role': 'system',
                 'content': 'You are a plant systems biologist. '
                            'Return ONLY valid JSON — no markdown, no preamble.'},
                {'role': 'user', 'content': prompt},
            ],
            'max_tokens': 2000,
            'temperature': 0.3,
            'response_format': {'type': 'json_object'},
        }
        if thinking:
            if _is_openai_provider(api_base):
                payload['reasoning'] = {'effort': reasoning_effort}
            else:
                payload['thinking'] = {'type': 'enabled'}
            payload['max_tokens'] = 8000  # thinking/reasoning + answer share budget
            payload.pop('response_format', None)
        label = f"{met_name} ({i}/{total})"
        print(f"  [{label}] {len(big_mods)} modules → analyzing …")
        try:
            body = _call_api(url, headers, payload, verbose=True, label=label)
            raw = body['choices'][0]['message']['content'].strip()
            parsed = _parse_json_response(raw, label)
            if parsed:
                if isinstance(parsed, list):
                    parsed = {'modules': parsed,
                              'metabolite': met_name,
                              'network_summary': ''}
                results[met_name] = parsed
                print(f"  [{label}]   ✓ {len(parsed.get('modules',[]))} modules "
                      f"characterized")
            else:
                print(f"  [{label}]   ⚠ parse failed, skipping")
        except Exception as e:
            print(f"  [{label}]   ✗ API error: {e}")
        if batch_delay > 0:
            time.sleep(batch_delay)
    return results


def _export_cross_metabolite_summary(isoform_presence, all_networks, out_prefix):
    """Export cross-metabolite summary: pleiotropic isoforms."""
    # Count metabolites per isoform
    pleio = [(iso, mets) for iso, mets in isoform_presence.items()
             if len(mets) >= 2]
    pleio.sort(key=lambda x: -len(x[1]))

    if not pleio:
        return

    path = f"{out_prefix}_cross_metabolite_isoforms.tsv"
    rows = []
    for iso, mets in pleio:
        rows.append({
            'Isoform': iso,
            'N_metabolites': len(mets),
            'Metabolites': ', '.join(sorted(mets)),
        })
    pd.DataFrame(rows).to_csv(path, sep='\t', index=False)
    print(f"\n[CROSS-METAB] → {path}  ({len(rows)} pleiotropic isoforms "
          f"appearing in ≥2 metabolite networks)")

    # Top-20 to console
    print(f"[CROSS-METAB]   top 15 pleiotropic isoforms:")
    for iso, mets in pleio[:15]:
        print(f"    {iso:22s}  {len(mets):2d} metabolites: "
              f"{', '.join(mets[:5])}{' …' if len(mets) > 5 else ''}")


def _plot_top_metabolite_networks(all_networks, all_ai, args):
    """Plot the top metabolite networks (most edges)."""
    # Rank metabolites by edge count
    ranked = []
    for met_name, net in all_networks.items():
        ranked.append((len(net['edges']), met_name))
    ranked.sort(reverse=True)

    # Plot top 12 (or fewer)
    top_n = min(12, len(ranked))
    print(f"\n[PLOT] Rendering top {top_n} metabolite networks …")

    for _, met_name in ranked[:top_n]:
        net = all_networks[met_name]
        ai = all_ai.get(met_name)
        safe_name = re.sub(r'[^a-zA-Z0-9_-]', '_', met_name)
        met_prefix = f"{args.out_prefix}_{safe_name}"

        # Export just this metabolite's files if not already done
        # (they are already exported, so R just needs to find them)
        modules = net['modules']
        has_mod = any(len(m) >= 2 for m in modules)
        if not has_mod:
            continue

        _run_r_plot(met_prefix, args.plot_script)


def _run_r_plot(out_prefix, plot_script=None):
    """Invoke Rscript to render the co-function network plot.

    Parameters
    ----------
    out_prefix : str
        Used to build input file paths and output prefix.
    plot_script : str or None
        Path to plot_cofunction_network.R.  If None, looks for
        ``plot_cofunction_network.R`` alongside this Python script.
    """
    if plot_script is None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        plot_script = os.path.join(script_dir, 'plot_cofunction_network.R')

    if not os.path.exists(plot_script):
        print(f"[PLOT] ERROR: R plot script not found: {plot_script}")
        print(f"[PLOT] Specify --plot-script PATH or ensure "
              f"plot_cofunction_network.R is present.")
        return

    nodes_path = f"{out_prefix}_cofunction_nodes.tsv"
    edges_path = f"{out_prefix}_cofunction_edges.tsv"
    modules_path = f"{out_prefix}_cofunction_modules.tsv"
    ai_json_path = f"{out_prefix}_network_ai_analysis.json"

    cmd = ['Rscript', plot_script, nodes_path, edges_path, modules_path]
    # Pass AI analysis JSON if it exists (5th positional arg)
    if os.path.exists(ai_json_path):
        cmd.append(ai_json_path)
        cmd.append(out_prefix)
    else:
        cmd.append('NONE')
        cmd.append(out_prefix)

    print(f"\n[PLOT] Rendering co-function network with R …")
    import subprocess
    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
    except FileNotFoundError:
        print("[PLOT] WARNING: Rscript not found — skipping network plots. "
              "(network tables/TSVs were still written).", file=sys.stderr)
        return
    if result.stdout:
        for line in result.stdout.strip().split('\n'):
            print(f"  {line}")
    if result.stderr:
        for line in result.stderr.strip().split('\n'):
            if 'mbcsToSbcs' not in line and 'conversion failure' not in line:
                print(f"  [R] {line}", file=sys.stderr)
    if result.returncode != 0:
        print(f"[PLOT] R exited with code {result.returncode}", file=sys.stderr)
    else:
        print(f"[PLOT] Network plots generated successfully")


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def write_output(all_rows, out_prefix):
    """Write the consolidated annotated table and summary statistics.

    Parameters
    ----------
    all_rows : list of dict
        All annotated rows across all metabolites.
    out_prefix : str
        Output file prefix.
    """
    if not all_rows:
        print("[OUTPUT] WARNING: No significant hits found — no output written.")
        return

    df = pd.DataFrame(all_rows)

    # Column order for TSV
    out_cols = [
        'Metabolites', 'Chr', 'Gene_start', 'Gene_end',
        'ref_Gene', 'Isoform', 'Pvalue',
        'Annotation', 'GO_description', 'PFAM',
    ]
    # Add AI columns if present
    if 'AI_link' in df.columns:
        out_cols.append('AI_link')
    if 'AI_confidence' in df.columns:
        out_cols.append('AI_confidence')

    avail_cols = [c for c in out_cols if c in df.columns]
    tsv_path = f"{out_prefix}_significant_isoforms.tsv"
    df[avail_cols].to_csv(tsv_path, sep='\t', index=False)
    print(f"\n[OUTPUT] → {tsv_path}  ({len(df)} rows × {len(avail_cols)} cols)")

    # --- Summary statistics ---
    stats_path = f"{out_prefix}_summary_stats.txt"

    total_hits = len(df)
    n_metabolites = df['Metabolites'].nunique()

    # TCONS vs LOC
    n_tcons = len(df[df['Isoform'].str.startswith('TCONS')])
    n_loc = len(df[df['Isoform'].str.startswith('LOC')])

    # Annotation coverage
    n_annotated = len(df[df['Annotation'].notna() & (df['Annotation'] != '')])
    n_go_desc = len(df[df['GO_description'].notna() & (df['GO_description'] != '')])
    n_pfam = len(df[df['PFAM'].notna() & (df['PFAM'] != '')])

    with open(stats_path, 'w') as f:
        f.write("# GWAS Significant Isoforms — Summary Statistics\n")
        f.write(f"# Generated by ISONET_LAYER1.py\n\n")
        f.write(f"Total significant isoforms:    {total_hits:,}\n")
        f.write(f"Unique metabolites:            {n_metabolites}\n")
        f.write(f"\n")
        f.write(f"Isoform type:\n")
        f.write(f"  TCONS (transcript):          {n_tcons:,}  ({100*n_tcons/max(total_hits,1):.1f}%)\n")
        f.write(f"  LOC (gene-level):            {n_loc:,}  ({100*n_loc/max(total_hits,1):.1f}%)\n")
        f.write(f"\n")
        f.write(f"Annotation coverage:\n")
        f.write(f"  With Description:            {n_annotated:,}  ({100*n_annotated/max(total_hits,1):.1f}%)\n")
        f.write(f"  With GO descriptions:        {n_go_desc:,}  ({100*n_go_desc/max(total_hits,1):.1f}%)\n")
        f.write(f"  With PFAM domains:           {n_pfam:,}  ({100*n_pfam/max(total_hits,1):.1f}%)\n")
        f.write(f"\n")
        f.write(f"Per-metabolite breakdown:\n")
        for metab, grp in df.groupby('Metabolites'):
            f.write(f"  {metab:40s}  {len(grp):6,} hits\n")

        # Top-20 most significant overall
        f.write(f"\nTop 20 most significant hits (overall):\n")
        top20 = df.sort_values('Pvalue').head(20)
        for _, row in top20.iterrows():
            f.write(f"  {row['Metabolites']:35s}  "
                    f"{row['Isoform']:20s}  "
                    f"p={row['Pvalue']:.2e}  "
                    f"{str(row.get('Annotation',''))[:70]}\n")

    print(f"[OUTPUT] → {stats_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="ISONET_LAYER1 — multi-metabolite GWAS to annotated isoform table + co-function network",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # --- Input sources (mutually exclusive modes for specifying QQMan files) ---
    qqman_group = parser.add_mutually_exclusive_group(required=True)
    qqman_group.add_argument(
        '--qqman-dir', dest='qqman_dir',
        help='Directory containing QQMan files (use with --qqman-pattern)',
    )
    qqman_group.add_argument(
        '--qqman-files', dest='qqman_files', nargs='+',
        help='Explicit list of QQMan file paths',
    )

    parser.add_argument(
        '--qqman-pattern', dest='qqman_pattern', default='*.pheno.qqman',
        help='Glob pattern to match QQMan files in --qqman-dir',
    )
    parser.add_argument(
        '--gtf', required=True,
        help='Path to StringTie combined GTF (TCONS→LOC mapping + coordinates)',
    )
    parser.add_argument(
        '--gff', default=None,
        help='Path to RefSeq GFF3 (gene coordinates for LOC entries)',
    )
    parser.add_argument(
        '--annotations', required=True,
        help='Path to eggNOG-mapper annotations file (TCONS + LOC query IDs)',
    )
    parser.add_argument(
        '--go-tab', dest='go_tab', default=None,
        help='Path to go.tab (GO ID → human-readable description)',
    )

    # --- Filtering ---
    parser.add_argument(
        '--pvalue-threshold', dest='pvalue_threshold', type=float, default=1e-5,
        help='Maximum p-value for a hit to be included in the output',
    )

    # --- Output ---
    parser.add_argument(
        '--out-prefix', dest='out_prefix', default='metabolites_significant',
        help='Prefix for all output files',
    )

    # --- AI enrichment (optional) ---
    parser.add_argument(
        '--api-key', dest='api_key', default=None,
        help='API key for AI enrichment (DeepSeek or OpenAI-compatible)',
    )
    parser.add_argument(
        '--api-base', dest='api_base', default='https://api.deepseek.com/v1',
        help='OpenAI-compatible API base URL (default: DeepSeek). '
             'Use https://api.openai.com/v1 for OpenAI models.',
    )
    parser.add_argument(
        '--model', dest='model', default='deepseek-v4-pro',
        help='Model name for AI enrichment (default: deepseek-v4-pro). '
             'Also: deepseek-chat, gpt-5.6 (OpenAI), gpt-5.5. '
             'For OpenAI reasoning depth use --thinking with --reasoning-effort.',
    )
    parser.add_argument(
        '--thinking', dest='thinking', action='store_true', default=False,
        help='Enable model reasoning before responding. '
             'For OpenAI models → reasoning.effort (controlled by --reasoning-effort). '
             'For DeepSeek models → thinking {type: enabled}. '
             'Improves hypothesis quality for complex regulatory analysis.',
    )
    parser.add_argument(
        '--reasoning-effort', dest='reasoning_effort', default='medium',
        choices=['low', 'medium', 'high'],
        help='OpenAI reasoning effort level when --thinking is enabled '
             '(default: medium). "low" = faster/cheaper, "high" = most thorough. '
             'Ignored for DeepSeek models.',
    )
    parser.add_argument(
        '--ai-system', dest='ai_system', default=None,
        help='Custom system instruction for the AI model '
             '(overrides the default plant molecular biologist persona)',
    )
    parser.add_argument(
        '--ai-rate-limit', dest='ai_rate_limit', type=float, default=0.1,
        help='Seconds between AI API calls (rate limiting)',
    )

    # --- Co-function network + AI regulatory analysis (optional) ---
    parser.add_argument(
        '--build-network', dest='build_network', action='store_true', default=False,
        help='Build co-function network from significant isoforms and '
             '(with --api-key) use AI to identify regulatory modules',
    )
    parser.add_argument(
        '--network-top-n', dest='network_top_n', type=int, default=300,
        help='Number of top significant isoforms to include in the network '
             '(taken from the pooled per-metabolite top-N set)',
    )
    parser.add_argument(
        '--network-per-metabolite-top', dest='network_per_metab_top',
        type=int, default=300,
        help='For each metabolite, take the top-N most significant isoforms '
             'from its QQMan file (no p-value cutoff).  Pooled across metabolites, '
             'deduplicated, then the global top --network-top-n are retained for '
             'the co-function network.',
    )
    parser.add_argument(
        '--network-min-shared', dest='network_min_shared', type=int, default=1,
        help='Minimum number of shared GO/PFAM terms to create an edge '
             'between two isoforms',
    )
    parser.add_argument(
        '--network-min-module-size', dest='network_min_module', type=int, default=3,
        help='Minimum module size for AI analysis (smaller modules reported but '
             'not sent to AI)',
    )
    parser.add_argument(
        '--plot-network', dest='plot_network', action='store_true', default=False,
        help='Generate publication-quality co-function network plot via R '
             '(requires Rscript with ggplot2, igraph, ggrepel)',
    )

    # --- Legacy network mode (preserved for backward compatibility) ---
    parser.add_argument(
        '--mode', choices=['table', 'network'], default='table',
        help='Output mode: "table" (default) for annotated isoform table, '
             '"network" for the legacy Cytoscape network mode',
    )
    # Legacy network arguments
    parser.add_argument(
        '--qqman', default=None,
        help='[legacy network mode] Path to a single QQMan file',
    )
    parser.add_argument(
        '--phenotype', default='Trait',
        help='[legacy network mode] Phenotype name (network center node)',
    )
    parser.add_argument(
        '--top-n', dest='top_n', type=int, default=300,
        help='[legacy network mode] Number of top significant unique genes to retain',
    )
    parser.add_argument(
        '--plot', action='store_true', default=False,
        help='[legacy network mode] Generate R network plots',
    )
    parser.add_argument(
        '--plot-script', dest='plot_script', default=None,
        help='Path to R plotting script (plot_network.R for --mode network, '
             'or plot_cofunction_network.R for --plot-network). '
             'Default: looks alongside this script.',
    )

    args = parser.parse_args()

    # ------------------------------------------------------------------
    # Legacy network mode
    # ------------------------------------------------------------------
    if args.mode == 'network':
        if not args.qqman:
            parser.error("--qqman is required for --mode network")
        run_legacy_network_mode(args)
        return

    # ------------------------------------------------------------------
    # Table mode (default) — multi-metabolite batch pipeline
    # ------------------------------------------------------------------
    run_table_mode(args)


def run_table_mode(args):
    """Multi-metabolite table-generation pipeline."""
    # --- Collect QQMan files ---
    if args.qqman_files:
        qqman_files = sorted(args.qqman_files)
    else:
        import glob
        pattern = os.path.join(args.qqman_dir, args.qqman_pattern)
        qqman_files = sorted(glob.glob(pattern))
        if not qqman_files:
            print(f"ERROR: No files matching '{pattern}'", file=sys.stderr)
            sys.exit(1)

    # Derive metabolite names from filenames
    # e.g. "2-methyl-malate.pheno.qqman" → "2-methyl-malate"
    metabolite_names = []
    for f in qqman_files:
        basename = os.path.basename(f)
        # Strip common suffixes
        name = basename
        for suffix in ['.pheno.qqman', '.qqman', '.pheno', '.txt']:
            if name.endswith(suffix):
                name = name[:-len(suffix)]
                break
        metabolite_names.append(name)

    print(f"[INPUT] {len(qqman_files)} QQMan files loaded")
    for i, (f, m) in enumerate(zip(qqman_files, metabolite_names)):
        print(f"  {i+1:3d}. {m:45s} ← {os.path.basename(f)}")

    # --- Load reference data ---
    t0 = time.time()

    tcons_map = parse_gtf_transcripts(args.gtf) if args.gtf else {}
    gff_gene_map = {}
    chr_map = {}
    if args.gff and os.path.exists(args.gff):
        gff_gene_map, chr_map = parse_gff(args.gff)

    emapper = load_emapper(args.annotations)

    go_map = {}
    if args.go_tab and os.path.exists(args.go_tab):
        go_map = load_go_map(args.go_tab)
    elif args.go_tab:
        print(f"[GO] WARNING: go.tab not found: {args.go_tab}")

    t1 = time.time()
    print(f"[LOAD] Reference data loaded in {t1 - t0:.1f}s\n")

    # --- Process each metabolite ---
    all_rows = []
    for qqman_path, metab_name in zip(qqman_files, metabolite_names):
        rows = process_qqman(
            qqman_path=qqman_path,
            metabolite_name=metab_name,
            pvalue_threshold=args.pvalue_threshold,
            tcons_map=tcons_map,
            gff_gene_map=gff_gene_map,
            emapper=emapper,
            go_map=go_map,
        )
        all_rows.extend(rows)

    total_hits = len(all_rows)
    print(f"\n[RESULT] {total_hits:,} significant hits across "
          f"{len(set(r['Metabolites'] for r in all_rows))} metabolites")
    if total_hits == 0:
        print("[RESULT] No hits to output. Try relaxing --pvalue-threshold.")
        return

    # --- Optional AI enrichment ---
    if args.api_key:
        all_rows = enrich_with_ai(
            all_rows,
            api_key=args.api_key,
            api_base=args.api_base,
            model=args.model,
            system_instruction=args.ai_system,
            batch_delay=args.ai_rate_limit,
            thinking=args.thinking,
            reasoning_effort=args.reasoning_effort,
        )

    # --- Write output ---
    write_output(all_rows, args.out_prefix)

    # --- Optional co-function network + AI regulatory analysis ---
    if args.build_network:
        run_network_pipeline(
            qqman_files, metabolite_names, args,
            tcons_map, gff_gene_map, emapper, go_map,
        )

    t2 = time.time()
    print(f"\n[DONE] Completed in {t2 - t0:.1f}s total")


# ---------------------------------------------------------------------------
# Legacy network mode (preserved from original IsoNet pipeline)
# ---------------------------------------------------------------------------

# COG categories kept for legacy network mode
COG_CATEGORIES = {
    'A': 'RNA processing and modification',
    'B': 'Chromatin structure and dynamics',
    'C': 'Energy production and conversion',
    'D': 'Cell cycle control, cell division, chromosome partitioning',
    'E': 'Amino acid transport and metabolism',
    'F': 'Nucleotide transport and metabolism',
    'G': 'Carbohydrate transport and metabolism',
    'H': 'Coenzyme transport and metabolism',
    'I': 'Lipid transport and metabolism',
    'J': 'Translation, ribosomal structure and biogenesis',
    'K': 'Transcription',
    'L': 'Replication, recombination and repair',
    'M': 'Cell wall/membrane/envelope biogenesis',
    'N': 'Cell motility',
    'O': 'Posttranslational modification, protein turnover, chaperones',
    'P': 'Inorganic ion transport and metabolism',
    'Q': 'Secondary metabolites biosynthesis, transport and catabolism',
    'R': 'General function prediction only',
    'S': 'Function unknown',
    'T': 'Signal transduction mechanisms',
    'U': 'Intracellular trafficking, secretion, and vesicular transport',
    'V': 'Defense mechanisms',
    'W': 'Extracellular structures',
    'X': 'Mobilome: prophages, transposons',
    'Y': 'Nuclear structure',
    'Z': 'Cytoskeleton',
}


def run_legacy_network_mode(args):
    """Original single-phenotype network-building pipeline.

    Preserved for backward compatibility — invokes the same logic as
    the original IsoNet pipeline but routed through --mode network.
    """
    from collections import Counter

    def _expand_cog(code):
        if not code or code == '-':
            return ''
        expanded = []
        for ch in code:
            name = COG_CATEGORIES.get(ch)
            if name:
                expanded.append(f"[{ch}] {name}")
        return '; '.join(expanded) if expanded else code

    def _safe_split(text, sep=','):
        if not text or text == '-' or pd.isna(text):
            return []
        return [t.strip() for t in str(text).split(sep) if t.strip() and t.strip() != '-']

    print("=== Legacy Network Mode ===")

    # 1. Parse QQMan
    print(f"[1/5] Parsing QQMan file: {args.qqman}")
    df = pd.read_csv(
        args.qqman, sep='\t', header=None,
        names=['original_id', 'chr', 'pos', 'pvalue'],
        dtype={'original_id': str, 'chr': str, 'pos': int, 'pvalue': float}
    )
    df['gene_id'] = df['original_id'].apply(extract_gene_id)
    df_dedup = df.loc[df.groupby('gene_id')['pvalue'].idxmin()].copy()
    df_dedup['neg_log10_p'] = -df_dedup['pvalue'].apply(lambda p: log10(max(p, 1e-300)))
    df_top = df_dedup.sort_values('pvalue').head(args.top_n).reset_index(drop=True)
    print(f"      Top {args.top_n} retained from {len(df_dedup):,} unique IDs")

    # 2. Load annotations
    print(f"[2/5] Loading annotations …")
    emapper = load_emapper(args.annotations)
    gff_fallback = {}
    if args.gff and os.path.exists(args.gff):
        gff_fallback, _ = parse_gff(args.gff)

    # 3. Annotate hits
    print(f"[3/5] Annotating {len(df_top)} hits …")
    for col in ['source', 'description', 'COG_category', 'COG_full_name',
                'Preferred_name', 'GOs', 'EC', 'KEGG_ko', 'KEGG_Pathway', 'PFAMs']:
        df_top[col] = ''

    matched_e_tcons, matched_e_loc, matched_gff, unmatched = 0, 0, 0, 0
    for idx, row in df_top.iterrows():
        gid = row['gene_id']
        if gid.startswith('TCONS'):
            df_top.at[idx, 'source'] = 'tcons'
            rec = emapper.get(gid)
            if rec:
                matched_e_tcons += 1
            else:
                unmatched += 1
                continue
        elif gid.startswith('LOC'):
            df_top.at[idx, 'source'] = 'loc'
            if gid in emapper:
                rec = emapper[gid]
                matched_e_loc += 1
            elif gff_fallback and gid in gff_fallback:
                g = gff_fallback[gid]
                df_top.at[idx, 'description'] = g['description']
                matched_gff += 1
                continue
            else:
                unmatched += 1
                continue
        else:
            unmatched += 1
            continue

        df_top.at[idx, 'description'] = rec.get('Description', '')
        df_top.at[idx, 'Preferred_name'] = rec.get('Preferred_name', '')
        df_top.at[idx, 'GOs'] = rec.get('GOs', '')
        df_top.at[idx, 'EC'] = rec.get('EC', '')
        df_top.at[idx, 'KEGG_ko'] = rec.get('KEGG_ko', '')
        df_top.at[idx, 'KEGG_Pathway'] = rec.get('KEGG_Pathway', '')
        df_top.at[idx, 'PFAMs'] = rec.get('PFAMs', '')
        cog = rec.get('COG_category', '')
        df_top.at[idx, 'COG_category'] = cog
        df_top.at[idx, 'COG_full_name'] = _expand_cog(cog)

    print(f"      Emapper TCONS: {matched_e_tcons}  |  Emapper LOC: {matched_e_loc}"
          f"  |  GFF: {matched_gff}  |  Unmatched: {unmatched}")

    # 4. Build network (simplified legacy output)
    print(f"[4/5] Building network …")
    nodes, edges = [], []
    seen_func = set()
    nodes.append({'id': args.phenotype, 'type': 'phenotype',
                  'display_name': args.phenotype.replace('_', ' ')})
    for _, row in df_top.iterrows():
        gid = row['gene_id']
        nodes.append({
            'id': gid, 'type': 'transcript' if row['source'] == 'tcons' else 'gene',
            'display_name': gid, 'description': row['description'],
            'chr': row['chr'], 'pos': row['pos'],
            'pvalue': f"{row['pvalue']:.4e}", 'neg_log10_p': row['neg_log10_p'],
        })
        edges.append({'source': args.phenotype, 'target': gid,
                      'weight': row['neg_log10_p'], 'type': 'association'})
        # Description as function term
        desc = row.get('description', '')
        if desc and desc != '-':
            fn = f"DESC:{desc[:120]}"
            if fn not in seen_func:
                seen_func.add(fn)
                nodes.append({'id': fn, 'type': 'description', 'display_name': desc[:120]})
            edges.append({'source': gid, 'target': fn, 'weight': 1.0, 'type': 'annotation'})

    nodes_df = pd.DataFrame(nodes)
    edges_df = pd.DataFrame(edges)
    print(f"      Nodes: {len(nodes_df)}  |  Edges: {len(edges_df)}")

    # 5. Export
    print(f"[5/5] Exporting …")
    nodes_df.to_csv(f"{args.out_prefix}_nodes.tsv", sep='\t', index=False)
    edges_df.to_csv(f"{args.out_prefix}_edges.tsv", sep='\t', index=False)
    df_top.to_csv(f"{args.out_prefix}_top_annotated.tsv", sep='\t', index=True, index_label='rank')
    print(f"      → {args.out_prefix}_nodes.tsv, {args.out_prefix}_edges.tsv")
    print("      Done.")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    main()
