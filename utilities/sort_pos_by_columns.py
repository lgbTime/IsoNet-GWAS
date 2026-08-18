import sys
import pandas as pd
def sort(args, ascend):
    df = pd.read_csv(args.df,header=None,sep="\t",index_col=None)
    if not args.col2:
        df.sort_values(by=[args.col1], inplace=True, ascending = ascend)
        df.to_csv("sorted_"+args.df, sep="\t",index=None,header=None)
    if args.col2 and args.col1:
        df.sort_values(by=[args.col1, args.col2], inplace=True, ascending = ascend)
        df.to_csv("sorted_"+args.df, sep="\t",index=None,header=None)

if __name__ == "__main__":
    import pandas as pd
    import argparse
    parser = argparse.ArgumentParser(description='sort file by col1 or (col1,col2) value')
    parser.add_argument('-f', dest='df', metavar='FILE', required=True,help='Path to first fasta')
    parser.add_argument('-c1', dest='col1', required=True, type=int, help='base on this columns to sort the file first')
    parser.add_argument('-c2', dest='col2', type=int, help='base on this columns to sort the file after col1')
    parser.add_argument('-t', dest='ascending', default="T", help='ascending setting T for True, F for False')
    args = parser.parse_args()
    if args.ascending == "T":
        ascend = True
    else:
        ascend = False
    sort(args, ascend)


