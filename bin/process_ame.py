#!/usr/bin/env python
import argparse
import os
import pandas as pd

def main(arg):
    annot = pd.read_csv(arg.annot, sep='\t', header=0, index_col=0)
    get_tf = annot.loc[:, ['Transcription factor']].drop_duplicates().to_dict(
        orient='dict')['Transcription factor']

    try:
        ame = pd.read_csv(arg.ame, sep='\t', header=0, comment='#')
    except:
        os.system("mv "+arg.ame+" "+arg.output)
        return
    ame['tf'] = [get_tf[m] for m in ame.motif_ID.values]
    ame.drop(['motif_DB', 'motif_alt_ID'], axis=1, inplace=True) 

    ame.to_csv(arg.output)
    return



if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='AME tsv to csv')

    parser.add_argument(
        '-f',
        '--ame',
        type=str,
        default='ame.tsv',
        help='Input file. (An AME output TSV) (default: ame.tsv)')

    parser.add_argument(
        '-a',
        '--annot',
        type=str,
        default=None,
        help='Annotation file for HOCOMOCO motifs (default: None)')

    parser.add_argument('-o',
                        '--output',
                        type=str,
                        default='output.csv',
                        help='Output CSV file (default: output.csv)')

    args = parser.parse_args()

    main(args)