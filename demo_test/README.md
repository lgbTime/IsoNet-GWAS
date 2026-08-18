# GLK Demo — gwasor GWAS → Layer-1 → Layer-2 on real lettuce data

A small, fully runnable example of the IsoNet-GWAS pipeline using the **GLK**
(GOLDEN2-LIKE transcription factor) project data: 134 lettuce accessions
genotyped for isoform PAVs plus the causal 1-bp deletion in GLK
(chr4:108,310,413), phenotyped for chlorophyll a+b.

Running the demo executes the whole flow with **no API key required**:

```
run_demo.sh
  ├── Step 1  gwasor.py pri_gwas      EMMAX GWAS (qq-norm → kinship → assoc → .qqman → plots)
  ├── Step 2  ISONET_LAYER1.py   significant isoforms → co-function network + modules
  └── Step 3  ISONET_LAYER2.py   expression-based expansion + regulatory ranking
```

## Quick start

```bash
# 1. point the tools in ../config.env at your PLINK/EMMAX/R/Python
#    (or export them, e.g. PLINK=/path/to/plink EMMAX=/path/to/emmax-intel64)

# 2. run the demo
bash run_demo.sh
```

## What you get

| Output | Content |
|--------|---------|
| `output/gwasout_all_chlorophy_a+b/` | EMMAX `.ps`/`.qqman` + Manhattan/QQ figures |
| `output/layer1/sig_chlorophyll_cofunction_{nodes,edges,modules}.tsv` | co-function network + modules (e.g. M1 protein-kinase module) |
| `output/layer1/*_network_*.pdf/png` | module-colored network figures (needs R) |
| `output/layer2/module_*_coexpression_edges.tsv` | co-expression edges of each expanded module |
| `output/layer2/module_*_regulatory_ranking.tsv` | regulators ranked by GWAS + expression score |
| `output/layer2/*_rising_hubs.tsv` | isoforms whose rank jumps after expansion |

## Options

- **AI module calling**: set `DEEPSEEK_API_KEY` (in `../config.env` or the
  environment) — Layer-1 then names modules and assigns high/medium/low
  confidence; Layer-2 can restrict expansion with `--ai-confidence high`.
- **Restrict expansion**: `DEMO_MODULES="M1 M3" bash run_demo.sh`.
- **Regenerate the demo data**: `bash build_demo_data.sh` (see below).

## The demo dataset (`data/`, 23 MB, committed on purpose)

The full GLK project files are far too large for GitHub (GTF 547 MB, GFF
223 MB, eggNOG 224 MB, expression matrix 1.5 GB), so `build_demo_data.sh`
built this subset from the original project:

| file | what it is |
|------|-----------|
| `addmode_Only4GLK.{bed,bim,fam,tfam}` | the original GLK genotype, kept **in full** (134 samples × ~63k isoform-PAV + SV markers) |
| `genetic_map.map` | tab-delimited isoform genetic map (chr, isoformID, cM, pos) |
| `all_chlorophy_a+b` | chlorophyll a+b phenotype |
| `go.tab` | GO-term description table (full) |
| `pop240GTF2ref.combined.gtf` | transcript lines only, for the p ≤ 1e-5 ∪ top-300 isoforms and all isoforms of their 127 parent genes |
| `Lsat_Salinas_v11_genomic.gff` | region lines (chromosome map) + 2,355 gene lines of the involved genes |
| `isoform_and_gene.emapper.annotations` | eggNOG rows of the involved isoforms/genes |
| `Pop_isoform_TPM_by_compareGTF.txt` | QC-filtered TPM matrix restricted to the involved isoforms + top-3000 most-expressed isoforms (3,289 × 240) |

Consequence: Layer-2 expansion candidates are limited to this subset
(≈3,300 isoforms), so the demo shows the mechanics on real data rather than
reproducing the paper-scale results.

## Notes

- The GWAS step needs EMMAX + PLINK (set in `../config.env`); R is used for
  phenotype qq-norm (bundled `../utilities/qqnorm_pheno.R`, base-R only)
  and for figures (qqman/GenABEL + ggplot2/igraph...).
- `gwasor.py` uses the `-map` argument (the tab-delimited isoform genetic
  map) to convert EMMAX `.ps` output into the `.qqman` format that carries
  each isoform's chromosome + physical position for the Manhattan plot.
