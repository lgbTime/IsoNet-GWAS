import pandas as pd
import argparse

def pav_xls2ped(pavxls,map_file):
    xls = pd.read_csv(pavxls, header=0, sep="\t", index_col=0)
    xls.drop(columns=["XLOC","refGene","classCode"], inplace=True)
    isomap = pd.read_csv(map_file, header=None, sep="\t", index_col=1)
    isoform = isomap.index

    xls_fil = xls.loc[isoform,]
    print(xls_fil)
    xls_fil = xls_fil.replace(1, "GG")
    xls_fil = xls_fil.replace(0, "AA").T
    insertName = xls_fil.index
    insert_zero = [0 for _ in range(len(insertName))]
    xls_fil.insert(0, "samples", insertName)
    xls_fil.insert(1, "zero1", insert_zero)
    xls_fil.insert(2, "zero2", insert_zero)
    xls_fil.insert(3, "zero3", insert_zero)
    xls_fil.insert(4, "zero4", insert_zero)
    xls_fil.to_csv(pavxls + ".ped", header=None, sep="\t", index=1)

def main():
    parser = argparse.ArgumentParser(description='Convert a PAV XLS file to PED format')
    parser.add_argument('pavxls', help='Path to the input PAV XLS file')
    parser.add_argument("map",    help= "genetics map use in gwas")
    args = parser.parse_args()
    pav_xls2ped(args.pavxls, args.map)

if __name__ == '__main__':
    main()
