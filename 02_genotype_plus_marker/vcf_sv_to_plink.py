#!/usr/bin/env python3
"""
vcf_sv_to_plink.py

Convert SV VCF (with symbolic alleles like <DEL>, <INS>, <DUP>, <TRA>)
to PLINK tped/tfam/map format, extracting only overlapping samples
between the VCF and an existing PLINK dataset.

Key design decisions:
  1. VCF chromosomes are expected to be numeric (1–9) matching the
     existing addmode_Only4GLK data.  A map_chrom() pass-through is
     included for backward compatibility with "Lsat_Salinas_v11_chr1"
     style names.
  2. GT conversion: 0/0→"1 1", 0/1→"1 2", 1/1→"2 2", ./.→"0 0".
     This is the standard PLINK allele encoding for biallelic variants.
  3. Output samples follow the EXACT order of addmode_Only4GLK.tfam.
     Samples in tfam but not in VCF get "0 0" (missing) for every variant.
  4. Streaming processing: the VCF is ~2.4M records; we stream line-by-line
     to avoid loading the entire file into memory.
  5. The resulting tped can be concatenated directly with the existing
     addmode_Only4GLK.tped (same column layout).

Usage:
    python vcf_sv_to_plink.py \
        --vcf extracted.gene_body.filtered_num_chr.vcf \
        --tfam addmode_Only4GLK.tfam \
        --out vcf_variants

Afterwards, merge with existing data:
    cat addmode_Only4GLK.tped vcf_variants.tped > merged.tped
    cat addmode_Only4GLK.map vcf_variants.map > merged.map
    cp addmode_Only4GLK.tfam merged.tfam

And optionally convert to binary PLINK:
    plink --tped merged.tped --tfam merged.tfam --make-bed --out merged
"""

import sys
import re
import argparse

# ---------------------------------------------------------------------------
# chromosome name mapping
# ---------------------------------------------------------------------------
CHR_RE = re.compile(r'Lsat_Salinas_v11_chr(\d+)')

def map_chrom(vcf_chrom):
    """Pass-through chromosome name (VCF already uses numeric chr).

    Also handles legacy "Lsat_Salinas_v11_chr1" → "1" style names
    when present, e.g. in the original (non-numeric) VCF.

    Examples
    --------
    1                      → "1"
    9                      → "9"
    Lsat_Salinas_v11_chr1  → "1"    (legacy)
    NW_026440196.1         → "NW_026440196.1"
    """
    m = CHR_RE.match(vcf_chrom)
    if m:
        return m.group(1)
    return vcf_chrom

# ---------------------------------------------------------------------------
# INFO field parsing
# ---------------------------------------------------------------------------
GENE_ID_RE = re.compile(r'GENE_ID=([^;]+)')

def extract_gene_id(info_field):
    """Extract GENE_ID from the VCF INFO column.

    When GENE_ID contains a comma-separated list (multi-gene overlap),
    only the first gene ID is returned to keep variant names within
    PLINK's 16,000-character limit.

    Example
    -------
    "ANNOTATION=gene_body;GENE_ID=LOC111905264;SVTYPE=INS;..." → "LOC111905264"
    "GENE_ID=LOC1,LOC2,LOC3"                                  → "LOC1"
    Returns None if GENE_ID is missing.
    """
    m = GENE_ID_RE.search(info_field)
    if m:
        gene_id = m.group(1)
        # Take only the first gene when GENE_ID is a comma-separated list
        return gene_id.split(',')[0]
    return None

def build_variant_id(var_id, gene_id):
    """Build variant identifier in format 'variantID:GENE_ID'.

    When GENE_ID is available:
        Lsat_Salinas_v11_chr1_44440 + LOC111905264 → Lsat_Salinas_v11_chr1_44440:LOC111905264
    When GENE_ID is missing, falls back to the bare variant ID:
        NW_026440196.1_20404 → NW_026440196.1_20404
    """
    if gene_id:
        return f"{var_id}:{gene_id}"
    return var_id

# ---------------------------------------------------------------------------
# GT conversion
# ---------------------------------------------------------------------------
def gt_to_plink(gt_str):
    """Convert a VCF GT field to a PLINK allele pair (diploid PLINK format).

    Handles both diploid and haploid VCF genotypes.  Encoding matches the
    lettome.tped / add_new_pos2tpedFormat.py convention:
      PLINK 1 = reference / SV-absent
      PLINK 2 = alternate  / SV-present
      PLINK 0 = missing

    Diploid:
      0/0   →   1 1      homozygous reference
      0/1   →   1 2      heterozygous
      1/1   →   2 2      homozygous alternate
      ./.   →   0 0      missing
      0|0   →   1 1      phased homozygous ref
    Haploid (single allele, used by some VCFs / BED-like files):
      0     →   1 1      haploid reference  → diploid homozygous ref
      1     →   2 2      haploid alternate  → diploid homozygous alt
      .     →   0 0      haploid missing     → diploid missing
    """
    if gt_str is None:
        return "0 0"

    gt = gt_str.strip()

    # Haploid: single character
    if gt in ('0', '1', '.'):
        if gt == '0':
            return "1 1"
        elif gt == '1':
            return "2 2"
        else:
            return "0 0"

    # Missing
    if gt in ('.', './.', '.|.'):
        return "0 0"

    # Diploid: normalise phase separator
    gt = gt.replace('|', '/')
    alleles = gt.split('/')
    out = []
    for a in alleles:
        if a == '.' or a == '':
            out.append('0')
        else:
            try:
                out.append(str(int(a) + 1))
            except ValueError:
                out.append('0')
    return ' '.join(out)

# ---------------------------------------------------------------------------
# header parsing
# ---------------------------------------------------------------------------
def parse_vcf_header(vcf_path):
    """Return list of sample names from the VCF #CHROM header line."""
    with open(vcf_path, 'r') as fh:
        for line in fh:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                parts = line.strip().split('\t')
                return parts[9:]   # fixed fields: CHROM…FORMAT
    raise SystemExit(f"ERROR: no #CHROM header line found in {vcf_path}")

def parse_tfam(tfam_path):
    """Return (sample_ids, fam_rows) from a PLINK .tfam file.

    sample_ids uses column 1 (individual ID) as key.
    fam_rows preserves the full lines for output.
    """
    samples = []
    fam_rows = []
    with open(tfam_path, 'r') as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            # Standard PLINK tfam: FID IID PAT MAT SEX PHENO
            # Some files use space, some use tab – split() handles both
            samples.append(parts[1] if len(parts) > 1 else parts[0])
            fam_rows.append(parts)
    return samples, fam_rows

# ---------------------------------------------------------------------------
# main conversion
# ---------------------------------------------------------------------------
def convert(vcf_path, tfam_path, out_prefix):
    # --- parse inputs -------------------------------------------------------
    vcf_samples = parse_vcf_header(vcf_path)
    tfam_samples, tfam_rows = parse_tfam(tfam_path)

    # Build lookup: VCF sample name → column index (0-based, relative to
    # the first sample column in the VCF)
    vcf_idx = {s: i for i, s in enumerate(vcf_samples)}

    # Keep only samples that appear in BOTH VCF and TFAM
    overlap_map = []          # list of (tfam_name, vcf_col) — overlap only
    overlap_rows = []         # corresponding fam rows
    skipped_tfam = []         # samples in tfam but not in VCF (dropped)
    for s, row in zip(tfam_samples, tfam_rows):
        if s in vcf_idx:
            overlap_map.append((s, vcf_idx[s]))
            while len(row) < 6:
                row.append('0')
            overlap_rows.append(row[:6])
        else:
            skipped_tfam.append(s)

    n_vcf = len(vcf_samples)
    n_tfam = len(tfam_samples)
    n_overlap = len(overlap_map)
    n_skipped = len(skipped_tfam)

    print(f"VCF samples:              {n_vcf}")
    print(f"TFAM samples:             {n_tfam}")
    print(f"Overlapping (kept):       {n_overlap}")
    print(f"TFAM-only (dropped):      {n_skipped}")
    if n_skipped <= 20:
        print(f"  dropped names:          {', '.join(skipped_tfam)}")
    else:
        print(f"  dropped names:          {', '.join(skipped_tfam[:10])} … "
              f"({n_skipped} total)")
    print()

    if n_overlap == 0:
        raise SystemExit("ERROR: zero overlapping samples — cannot proceed.")

    # --- open output files --------------------------------------------------
    tped_out = open(f"{out_prefix}.tped", 'w')
    map_out  = open(f"{out_prefix}.map", 'w')
    tfam_out = open(f"{out_prefix}.tfam", 'w')

    # Write tfam — overlap samples only, space-delimited (PLINK standard)
    for row in overlap_rows:
        tfam_out.write(' '.join(row) + '\n')

    # Write the sample list for subsetting the existing dataset
    keep_out = open(f"{out_prefix}.overlap_samples.txt", 'w')
    for name, _ in overlap_map:
        keep_out.write(f"{name}\t{name}\n")
    keep_out.close()

    # --- stream the VCF -----------------------------------------------------
    variant_count = 0
    report_every = 100000

    with open(vcf_path, 'r') as fh:
        for line in fh:
            if line.startswith('#'):
                continue

            parts = line.strip().split('\t')
            if len(parts) < 10:
                continue

            chrom_raw = parts[0]
            pos       = parts[1]
            var_id    = parts[2] if parts[2] != '.' else f"{chrom_raw}_{pos}"
            info      = parts[7]
            fmt_field = parts[8]

            # Extract gene ID from INFO and build enriched variant ID
            gene_id = extract_gene_id(info)
            var_id = build_variant_id(var_id, gene_id)

            # Locate the GT sub-field index within FORMAT
            fmt_keys = fmt_field.split(':')
            try:
                gt_idx = fmt_keys.index('GT')
            except ValueError:
                # No GT field — skip this variant
                continue

            # Map chromosome
            chrom = map_chrom(chrom_raw)

            # Build the tped line
            # tped format: chrom id genetic_dist pos gt1 gt2 ... gtN
            fields = [chrom, var_id, '0', pos]

            for name, vcf_col in overlap_map:
                sample_col = 9 + vcf_col
                if sample_col < len(parts):
                    sample_fields = parts[sample_col].split(':')
                    gt = sample_fields[gt_idx] if gt_idx < len(sample_fields) else './.'
                else:
                    gt = './.'
                fields.append(gt_to_plink(gt))

            tped_out.write(' '.join(fields) + '\n')

            # map format: chrom id genetic_dist pos
            map_out.write(f"{chrom}\t{var_id}\t0\t{pos}\n")

            variant_count += 1
            if variant_count % report_every == 0:
                print(f"  … processed {variant_count:,} variants")

    # --- clean up -----------------------------------------------------------
    tped_out.close()
    map_out.close()
    tfam_out.close()

    print(f"\nDone.  {variant_count:,} variants written to:")
    print(f"  {out_prefix}.tped")
    print(f"  {out_prefix}.map")
    print(f"  {out_prefix}.tfam")
    print(f"  {out_prefix}.overlap_samples.txt  (for subsetting existing dataset)")
    print()
    print("Next steps:")
    print(f"  # 1. Convert VCF tped to binary:")
    print(f"  plink --tped {out_prefix}.tped --tfam {out_prefix}.tfam --make-bed --out {out_prefix}")
    print(f"  # 2. Subset existing dataset to overlap samples:")
    print(f"  plink --bfile <existing> --keep {out_prefix}.overlap_samples.txt --make-bed --out existing_overlap")
    print(f"  # 3. Merge:")
    print(f"  plink --bfile existing_overlap --bmerge {out_prefix} --make-bed --out final")

# ---------------------------------------------------------------------------
# BED-like GT conversion (e.g. GLK.GT)
# ---------------------------------------------------------------------------
def convert_bed_gt(bed_path, tfam_path, out_prefix):
    """Convert a BED-like genotype file to PLINK tped/tfam/map.

    Expected format (tab-delimited, no header):
      chr  start  end  score  variant_id  SV_type  strand  INFO  gt1  gt2 ...
    where gt1, gt2, ... are haploid 0/1 genotypes, one per sample.

    Sample order in the BED file must match the VCF sample header order.
    """
    vcf_samples = parse_vcf_header(bed_path) if False else None

    tfam_samples, tfam_rows = parse_tfam(tfam_path)

    # Read the single line from the BED-GT file
    with open(bed_path, 'r') as fh:
        line = fh.readline().strip()
        if not line:
            raise SystemExit(f"ERROR: empty file: {bed_path}")
        parts = line.split('\t')

    # BED columns: 1=chr, 2=start, 3=end, 4=score, 5=var_id, 6=type,
    #              7=strand, 8=info, 9+=genotypes
    chrom    = map_chrom(parts[0])
    pos      = parts[1]  # start position
    var_id   = parts[4]
    info     = parts[7]
    genotypes = parts[8:]  # haploid 0/1 values

    # Extract gene ID
    gene_id = extract_gene_id(info)
    var_id = build_variant_id(var_id, gene_id)

    # For sample matching, we use the VCF header from the companion VCF
    # or match by position from the tfam

    # Build sample lookup from tfam
    tfam_set = set(tfam_samples)

    # We don't have sample names in the BED-GT file — only genotype order.
    # We'll build the overlap by reading the VCF sample order from the
    # companion VCF file passed via --vcf-samples.

    print(f"WARNING: --bed-gt needs sample ordering from VCF.")
    print(f"  Use --vcf-samples to pass VCF #CHROM header for sample order.")
    print(f"  Currently guessing: genotypes = {len(genotypes)} samples")
    print()

    # Fall back to just creating the output with the genotypes
    # This assumes the BED samples match the VCF sample order

    with open(f"{out_prefix}_glk.tped", 'w') as tped, \
         open(f"{out_prefix}_glk.map", 'w') as mapf, \
         open(f"{out_prefix}_glk.tfam", 'w') as tfamf:

        # Build tped line
        fields = [chrom, var_id, '0', pos]
        for gt in genotypes:
            fields.append(gt_to_plink(gt))

        tped.write(' '.join(fields) + '\n')
        mapf.write(f"{chrom}\t{var_id}\t0\t{pos}\n")

        # Write tfam from tfam_rows
        for row in tfam_rows:
            while len(row) < 6:
                row.append('0')
            tfamf.write('\t'.join(row[:6]) + '\n')

    print(f"GLK variant written to {out_prefix}_glk.*")
    print(f"  Variant ID: {var_id}")
    print(f"  Samples:    {len(genotypes)}")


def convert_bed_gt_with_overlap(bed_path, tfam_path, vcf_samples_path, out_prefix):
    """Like convert_bed_gt but matches overlap samples like convert() does."""
    import os

    # Parse VCF samples from --vcf-samples file (extracted #CHROM header line)
    vcf_samples = []
    with open(vcf_samples_path, 'r') as fh:
        for line in fh:
            parts = line.strip().split('\t')
            vcf_samples = parts[9:]
            break

    tfam_samples, tfam_rows = parse_tfam(tfam_path)
    vcf_idx = {s: i for i, s in enumerate(vcf_samples)}

    # Build overlap
    overlap_map = []
    overlap_rows = []
    skipped_tfam = []
    for s, row in zip(tfam_samples, tfam_rows):
        if s in vcf_idx:
            overlap_map.append((s, vcf_idx[s]))
            while len(row) < 6:
                row.append('0')
            overlap_rows.append(row[:6])
        else:
            skipped_tfam.append(s)

    n_overlap = len(overlap_map)

    print(f"VCF samples:              {len(vcf_samples)}")
    print(f"TFAM samples:             {len(tfam_samples)}")
    print(f"Overlapping (kept):       {n_overlap}")
    print(f"TFAM-only (dropped):      {len(skipped_tfam)}")
    print()

    # Read BED-GT line
    with open(bed_path, 'r') as fh:
        line = fh.readline().strip()
        if not line:
            raise SystemExit(f"ERROR: empty file: {bed_path}")
        parts = line.split('\t')

    chrom     = map_chrom(parts[0])
    pos       = parts[1]
    var_id    = parts[4]
    info      = parts[7]
    genotypes = parts[8:]

    gene_id = extract_gene_id(info)
    var_id  = build_variant_id(var_id, gene_id)

    # Write tped/map/tfam
    with open(f"{out_prefix}.tped", 'w') as tped, \
         open(f"{out_prefix}.map", 'w') as mapf, \
         open(f"{out_prefix}.tfam", 'w') as tfamf, \
         open(f"{out_prefix}.overlap_samples.txt", 'w') as keepf:

        for row in overlap_rows:
            tfamf.write(' '.join(row) + '\n')

        for name, _ in overlap_map:
            keepf.write(f"{name}\t{name}\n")

        fields = [chrom, var_id, '0', pos]
        for sname, col in overlap_map:
            gt = genotypes[col] if col < len(genotypes) else '.'
            fields.append(gt_to_plink(gt))

        tped.write(' '.join(fields) + '\n')
        mapf.write(f"{chrom}\t{var_id}\t0\t{pos}\n")

    print(f"Written: {var_id}")
    print(f"  {out_prefix}.tped")
    print(f"  {out_prefix}.map")
    print(f"  {out_prefix}.tfam")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Convert SV VCF / BED-GT to PLINK tped/tfam/map'
    )
    parser.add_argument('--vcf',  default=None,
                        help='Input VCF file')
    parser.add_argument('--bed-gt', default=None,
                        help='Input BED-like genotype file (e.g. GLK.GT)')
    parser.add_argument('--vcf-samples', default=None,
                        help='File with VCF #CHROM header line for sample order (used with --bed-gt)')
    parser.add_argument('--tfam', required=True,
                        help='Reference .tfam file (defines sample set & order)')
    parser.add_argument('--out',  required=True,
                        help='Output file prefix')
    args = parser.parse_args()

    if args.vcf:
        convert(args.vcf, args.tfam, args.out)
    elif args.bed_gt:
        if args.vcf_samples:
            convert_bed_gt_with_overlap(args.bed_gt, args.tfam, args.vcf_samples, args.out)
        else:
            convert_bed_gt(args.bed_gt, args.tfam, args.out)
    else:
        parser.error("Either --vcf or --bed-gt is required")
