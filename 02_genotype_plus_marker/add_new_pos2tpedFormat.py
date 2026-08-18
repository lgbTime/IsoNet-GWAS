import pandas as pd
import argparse

# Create the argument parser
parser = argparse.ArgumentParser(description='Process genotypes for a specific gene.')
parser.add_argument('chrom', type=int, help='chromosome')
parser.add_argument('pos', type=int, help='position')
parser.add_argument('gene', type=str, help='gene name')
parser.add_argument('--tped', type=str, default='tped', help='path to TPED file')
parser.add_argument('--tfam', type=str, default='tfam', help='path to TFAM file')
parser.add_argument('--map', type=str, default='map', help='path to map file')

parser.add_argument('--gt', type=str, default='gt', help='path to genotype file, two columns (no head)')
parser.add_argument('--out', type=str, default='out', help='new tped name preffix')
# Parse the command-line arguments
args = parser.parse_args()

# Read the data files
tped = pd.read_csv(args.tped, header=None, sep=" ", index_col=None)
tfam = pd.read_csv(args.tfam, header=None, sep=" ", index_col=0)
mapdf = pd.read_csv(args.map, header=None, sep="\t", index_col=None)
gene_gts = pd.read_csv(args.gt, header=None, sep="\t", index_col=0)
gene_gts = gene_gts.loc[tfam.index, :]  # sort samples order
gts = list(set(gene_gts[1].tolist()))

def uniq_pos(tped,chrom, pos):
    while pos in tped[tped[0]==chrom][3].tolist():
        pos += 1
    return pos
    
if len(gts) < 3:
    gt_dict = dict(zip(gts, ["1 1", "2 2"]))
else:
    gt_dict = dict(zip(gts, ["1 1", "2 2", "1 2"]))

# Generate the output line
pos = uniq_pos(tped, args.chrom, args.pos)
out4map = f"{args.chrom} {args.gene} 0 {pos}"
out4tped =f"{args.chrom} {args.gene} 0 {pos}"
for r in gene_gts.index:
    gt = gene_gts.loc[r, 1]
    out4tped += " " + gt_dict[gt]
print(out4tped)
print(f"mapdf len:{len(mapdf)}")

## add new pos and genotype ##
linelst = out4tped.split(" ")
tped.loc[len(tped)] = linelst
tped.sort_values(by=[0,3],inplace=True) ## sort them by chrom, pos
tped.to_csv(args.out+".tped",header=None,sep=" ",index=None)

## add new gene map info ##
outmaplst = out4map.split(" ")
mapdf.loc[len(mapdf)] = outmaplst
mapdf.sort_values(by=[0,3]).to_csv(args.out+".map",header=None,sep=" ",index=None)


