# IGWAS — Isoform-PAV GWAS coupled with a two-layer network

**IGWAS** is a pangenome isoform-presence/absence-variation (PAV) GWAS framework
that couples genetic association with a **two-layer network** to move from
"significant isoform" to "regulatory module and hub isoform":

- **Layer 1 — co-function network**: significant isoforms (PAV + gene-body SV
  markers) are mapped to parent genes, annotated with eggNOG (GO / PFAM / KEGG /
  COG / descriptions), and wired into a co-function network in which two
  isoforms are connected when they share functional terms. Network modules are
  detected, and an optional LLM step (DeepSeek / OpenAI) names modules and
  assigns regulatory confidence (high / medium / low).
- **Layer 2 — expression network**: high-confidence Layer-1 modules are expanded
  with isoform-level co-expression (Spearman correlation + FDR), regulators are
  ranked by a combined GWAS + expression score, and "rising hubs" — isoforms
  whose rank jumps after expansion — are reported as candidate master
  regulators.

```
  isoform PAV matrix ─┐
  per-chr isoform BED ─┤ Step 1: unique-position genetic map ──► PLINK dataset
                      └► PAV -> PED -> bed/tped
                                        │
  gene-body SV VCF ───┐ Step 2: merge extra markers (optional)
  GLK PAV VCF ────────┘                    │
                                           ▼
  phenotype per trait ──► Step 3: EMMAX GWAS ──► *.ps ──► *.qqman (+ QQ/Manhattan)
                                           │
                                           ▼
                     Step 4: LAYER 1  co-function network (+ AI modules)
                                           │
                                           ▼
                     Step 5: LAYER 2  expression expansion ─► rising hubs
```

## Repository layout

```
IGWAS_GITHUB/
├── 00_master_pipeline.sh        # end-to-end runner (steps 1-5)
├── config.env                   # all tool paths, data locations, AI keys
├── 01_isoform_genetic_map/      # Step 1: PAV -> genetic map -> PLINK
│   ├── bed2map/                 #   unique-position isoform map
│   └── pav2ped/                 #   PAV matrix -> PED -> bed/tped
├── 02_genotype_plus_marker/     # Step 2: merge extra marker VCFs
├── 03_gwas/                     # Step 3: EMMAX GWAS + QQ/Manhattan
├── 04_layer1_cofunction_network/# Step 4: LAYER 1 co-function network + AI
├── 05_layer2_expression_network/# Step 5: LAYER 2 expression expansion
├── demo_test/                   # runnable demo: GLK data, gwasor→Layer1→Layer2
│   ├── run_demo.sh              #   one-command run (no API key required)
│   ├── build_demo_data.sh       #   rebuild the 23 MB demo subset from the full project
│   └── data/                    #   demo inputs (committed on purpose, see demo_test/README.md)
├── resources/                   # static annotation resource: go.tab (GO→description, from eggNOG-mapper)
└── utilities/                   # small helpers (incl. base-R phenotype qq-norm script)
```

## Requirements

**Required** (full pipeline: analysis + all figures):

| Required | Used in | Notes |
|----------|---------|-------|
| Python 3 + `pandas` | all steps | Layer 2 also needs `numpy`, `scipy`, `networkx` |
| R + Rscript | steps 3-5 | phenotype qq-norm, QQ/Manhattan and network figures |
| R packages: `qqman` `GenABEL` | step 3 | QQ/Manhattan figures + genomic-control lambda |
| R packages: `data.table` `ggplot2` `ggforce` `ggrepel` `igraph` `RColorBrewer` `viridisLite` `scales` | steps 4-5 | two-layer network node/module/co-expression figures |
| PLINK 1.9 | steps 1-3 | format conversion / merging |
| EMMAX (`emmax-intel64` + `emmax-kin-intel64`) | step 3 | association engine |

Phenotype qq-norm only needs base R (`qqnorm_pheno.R` is bundled in
`utilities/`); all other R packages are only used for figures — if you only
want the TSV analysis tables, the scripts skip plotting automatically (network
construction and `.qqman` output are unaffected).

**Optional**:

| Optional | Used in | Behaviour when missing |
|----------|---------|------------------------|
| flashpca2 | step 3 PCA covariates (`pca` subcommand) | not used by the default flow, can be ignored |
| DeepSeek API key | AI module calling in steps 4-5 | AI skipped, network construction unaffected |

All external tools are resolved from `config.env`; no hard-coded machine paths
are used anywhere in the repository.

## Quick start

1. **Prepare data** — create `data/` (git-ignored) with symlinks or copies of:
   - per-chromosome isoform BED files (col1=chr, col2=TSS, col3=TES, col4=isoform ID)
   - the isoform PAV matrix (`pav.xls`; rows=isoforms, sample columns, plus
     `XLOC`, `refGene`, `classCode` columns)
   - optional extra marker VCFs (gene-body SVs, GLK PAV)
   - phenotype files (one per trait; two columns: family ID, value)
   - `pop240GTF2ref.combined.gtf` (StringTie transcript GTF),
     `Lsat_Salinas_v11_genomic.gff` (reference annotation),
     `isoform_and_gene.emapper.annotations` — **eggNOG-mapper hit table, see
     "eggNOG annotation" below** (config `ANNOTATIONS`)
   - isoform TPM expression matrix (`Pop_isoform_TPM_by_compareGTF.txt`)

   `go.tab` (GO term description table) is **provided** in `resources/go.tab`
   and does not need to be prepared (config `GO_TAB` points to it by default).

   The isoform genetic map is **not** an input you prepare — Step 1 generates it
   automatically from the per-chromosome BED files (`$DATA_DIR/map/genetic_map.map`,
   one unique coordinate per isoform). It is then required by Step 3 for the
   `.qqman` conversion and Manhattan plotting, and by Step 1b for the PAV→PED
   conversion.
2. **Edit `config.env`** — set tool paths, data file locations, trait list and
   (optionally) `DEEPSEEK_API_KEY`.
3. **Run**:
   ```bash
   ./00_master_pipeline.sh          # steps 1-5
   ./00_master_pipeline.sh 3 5      # only GWAS + networks
   ```
   Each step can also be run individually from its own directory, e.g.
   `bash 04_layer1_cofunction_network/run_layer1.sh`.

**Want to see it run?** The repository ships a runnable demo (`demo_test/`,
23 MB of real GLK data, no API key required):

```bash
cd demo_test
bash run_demo.sh        # gwasor GWAS → Layer-1 → Layer-2, output in demo_test/output/
```

## Prerequisite — eggNOG annotation of isoforms and genes

Layer 1 and Layer 2 build the co-function network from functional annotations,
so **every isoform (TCONS) and its parent gene (LOC) must first be annotated
with [eggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper)** to produce
the hit table `isoform_and_gene.emapper.annotations` (config `ANNOTATIONS`).

**How to produce it:**

1. Merge the protein sequences of all isoforms and all genes into one FASTA
   (e.g. `hit_swiss.pep.fa` — isoform proteins + gene proteins).
2. Run eggNOG-mapper (this is the exact command used for the lettuce data,
   eggNOG-mapper 2.1.12):
   ```bash
   emapper.py -i hit_swiss.pep.fa -o hit_swiss.pep.fa --cpu 20
   # output: hit_swiss.pep.fa.emapper.annotations
   cp hit_swiss.pep.fa.emapper.annotations isoform_and_gene.emapper.annotations
   ```
3. Expected format: `##` comment lines, one tab-separated row per query, keyed
   by the query ID (`TCONS_*` for isoforms, `LOC*` for genes). The scripts look
   up annotations by either ID, so both must be present.

**`go.tab`** (GO ID → description table, used to render human-readable GO
terms on nodes/modules) comes from the eggNOG-mapper data directory and is
**provided in this repository** at `resources/go.tab` (config `GO_TAB` points
there by default; 3 columns: `GO`, `Description`, `level`).

## Step-by-step

### Step 1 — Isoform genetic map and PAV genotype (`01_isoform_genetic_map/`)

**The isoform genetic map is a core requirement of the whole pipeline** — every
isoform must have a unique chromosome coordinate because the Manhattan plot
(and the `.qqman` conversion in Step 3) needs the position of each isoform
marker. The map is produced here, then consumed by Step 3 (`-map` argument of
`gwasor.py`, config `MAPFILE`) and reused by the PLINK `.ped`/`.map` conversion.

- `bed2map/get_uniq_pos.sh` — deduplicates isoform TSS positions
  (`reassign_repeated_startpos2_v2.py`) so every isoform has a unique
  coordinate, then concatenates chromosomes into `genetic_map.bed` and a
  PLINK-ready `genetic_map.map`.

  `genetic_map.map` format (tab-separated, one isoform per line):

  | col | field | example |
  |-----|-------|---------|
  | 1 | chromosome | `Lsat_Salinas_v11_chr1` |
  | 2 | isoform ID | `TCONS_00004469` (or `LOC111905264`) |
  | 3 | genetic position (cM) | `0` |
  | 4 | physical position (bp) | `43133` |

  Only isoforms present in this map are kept in the PAV genotype (Step 1b) and
  plotted in the Manhattan figures.
- `pav2ped/plink.sh` — converts the isoform PAV matrix (0/1) to a PLINK PED
  (`pav_xls2ped.py`), then to binary PLINK and to transpose format
  (`--recode12 --transpose`), the genotype layout used by the GWAS step.

### Step 2 — Plus-marker genotype (`02_genotype_plus_marker/`)
- `vcf_sv_to_plink.py` — streams a structural-variant VCF (symbolic alleles,
  e.g. `<DEL>`, `<INS>`) into PLINK tped/tfam, keeping only samples shared
  with the PAV dataset; 0/0→"1 1", 0/1→"1 2", 1/1→"2 2".
- `run.sh` — merges the PAV dataset with the gene-body SV VCF and the GLK PAV
  VCF into one binary + text PLINK dataset (`lettome_iso_plus_gene_variants`).
- `add_new_pos2tpedFormat.py` — inserts one extra variant (e.g. a functional
  SNP genotype) at a given chromosome position into an existing tped/map.

### Step 3 — GWAS (`03_gwas/`)
- `gwasor.py` — EMMAX wrapper with subcommands:
  - `pri_gwas` : full chain — phenotype keep → quantile normalization →
    genotype subset (MAF 0.05) → kinship (BN/IBS) → association → `.ps` →
    `.qqman` → QQ/Manhattan plots
  - `fine_gwas`, `emmax-kin`, `emmax_genotype`, `pca`, `ps2plotdata`, `plot`

  **`-map` is mandatory for isoform-level GWAS**: `gwasor.py` needs the isoform
  genetic map (Step 1 output) to convert the EMMAX `.ps` result into the
  `.qqman` format (isoform ID + chromosome + physical position + p-value) —
  without the map there is nothing to plot. `gwas.sh` passes
  `$DATA_DIR/map/genetic_map.map` automatically (config `MAPFILE`).
- `gwas.sh` / `gwas_norm.sh` — loop over `$TRAITS` (e.g. pigment traits).
- `qqman4Pop_lecttuce.R` / `qqman4Pop_lecttuce_pdf.R` — QQ and Manhattan
  figures with genomic-control lambda; `snpsOfInterest` can be filled to
  highlight hits.

### Step 4 — Layer 1: co-function network (`04_layer1_cofunction_network/`)
- `IGWAS_DSNET_LAYER1.py` — significant isoforms (p < `P_VALUE_THRESHOLD`) from
  all `*.pheno.qqman` files → parent-gene mapping (GTF/GFF) → eggNOG annotation
  → co-function edges (shared GO/PFAM/KEGG/COG/description terms) → module
  detection → optional LLM module naming and confidence scoring
  (`--api-key`, requires `DEEPSEEK_API_KEY` in config.env).
- `build_v1_concentric.py` + `plot_layer1_network.R` — V1 concentric-circle
  network figure (full / co-function / with-orphans views).
- `plot_cofunction_network.R` — module-colored co-function plots, invoked
  automatically by `IGWAS_DSNET_LAYER1.py`.

### Step 5 — Layer 2: expression network (`05_layer2_expression_network/`)
- `prefilter_expression_matrix.py` — removes low-expression isoforms from the
  TPM matrix (mean TPM ≥ 0.5 **or** TPM ≥ 1 in ≥ 20 % of samples).
- `IGWAS_DSNET_LAYER2.py` — expands each AI-confirmed high-confidence module
  by expression co-variation, ranks regulators by a combined GWAS + expression
  score, and reports rising hubs.
- `plot_layer2_network.R` — per-module network / evidence-scatter / kME plots.

## Notes

- **Metabolome application**: the original project also applied Layer-1 to
  many metabolite GWAS files at once (`IGWAS_DSNET_V2.py`, a thin variant of
  `IGWAS_DSNET_LAYER1.py` with OpenAI/DeepSeek dual-provider AI support). It
  is not needed for the two-layer flow and is therefore not included here.
- **API keys**: the LLM step is optional. If `DEEPSEEK_API_KEY` is empty, the
  network construction and plotting still run; only the AI module-calling is
  skipped. Never commit a real key — `config.env` is designed to be filled in
  locally.
- **Data**: all inputs are referenced through `config.env` and are expected to
  live in the git-ignored `data/` directory; nothing in this repository depends
  on machine-specific paths.
- **Isoform PAV** is the core marker type; gene-body SVs and a GLK PAV variant
  are merged in Step 2 as examples of "plus markers" (see `run.sh` for the GT
  conversion table).

## License

TBD — please add a license before public release.
