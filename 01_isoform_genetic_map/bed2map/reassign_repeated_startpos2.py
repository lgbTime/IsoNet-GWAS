import pandas as pd
import argparse
def auto_rmrepeat_pos(df):
    """
    we need a uniq Isoform position file to generate Genetics map,
    here we use the each Isoform TSS as pos,but some Isoform may have the same position,
    then we try to use the TTS as pos if the TTS is unique.
    or we add offset to TSS to diff the pos since not unique TTS and TSS found.
    the new out_bed may still have few repeat pos, we should open the result file to manul correct them
    Last we must sort the file again by the chromosome and second columns
    """
    new_repeats = df[1].value_counts()
    for repeat_value, count_repeat in new_repeats.items():
        if count_repeat >= 2:
            tmpdf2 = df[df[1] == repeat_value]
            unique_positions = set(tmpdf2[2])
            
            if len(unique_positions) == tmpdf2.shape[0]:
                df.loc[tmpdf2.index, 1] = tmpdf2[2].tolist()
            
            else:
                new_positions = []
                for j, index in enumerate(tmpdf2.index):
                    offset = 1
                    while True:
                        new_position = int(df.loc[index, 1]) + j + offset
                        if new_position not in df[1].tolist() and new_position not in new_positions:
                            new_positions.append(new_position)
                            break
                        offset += 1
                df.loc[tmpdf2.index, 1] = new_positions
    return df


def main():
    parser = argparse.ArgumentParser(description ='To get a isoform map, we need a uniq Isoform position,here we use the TSS as position')
    parser.add_argument('bedfile', help='Path to the input bed file')
    args = parser.parse_args()

    bed = pd.read_csv(args.bedfile, sep="\t", index_col=None, header=None)
    bed.index = bed[3]
    for i in range(5):
        if len(set(bed[1].tolist())) != bed.shape[0]:
            bed = auto_rmrepeat_pos(bed)
    bed = bed.sort_values(by=[0, 1])
    bed.to_csv("out_" + args.bedfile, header=None, sep="\t", index=None)

if __name__ == '__main__':
    main()
