#!/usr/bin/env python
import argparse

import pandas as pd

def main(arg):
    annot = pd.read_csv(arg.annot, sep='\t', header=0, index_col=0)
    get_tf = annot.loc[:, ['Transcription factor']].drop_duplicates().to_dict(
        orient='dict')['Transcription factor']
    # annot.index = annot.index.str.replace('_MOUSE.*', '')
    # get_annot = annot.loc[:, ['Transcription factor']].drop_duplicates(
    # ).reset_index().set_index('Transcription factor').to_dict(
    #     orient='dict')['Model']

    fimo = pd.read_csv(arg.fimo, sep='\t', header=0, comment='#')
    fimo['tf'] = [get_tf[m] for m in fimo.motif_id.values]
    fimo = fimo.join(
        fimo['sequence_name'].str.split(
            '::', 1, expand=True).rename(columns={
                0: 'target_gene',
                1: 'location'
            }))
    fimo.drop(['sequence_name', 'motif_alt_id'], axis=1, inplace=True) 

    fimo.to_csv(arg.output, index=False)
    return



if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='FIMO tsv to xlsx')

    parser.add_argument(
        '-f',
        '--fimo',
        type=str,
        default='fimo.tsv',
        help='Input file. (A FIMO output TSV) (default: fimo.tsv)')

    parser.add_argument(
        '-a',
        '--annot',
        type=str,
        default=None,
        help='Annotation file for HOCOMOCO motifs (default: None)')

    parser.add_argument('-o',
                        '--output',
                        type=str,
                        default='output.xlsx',
                        help='Output Excel file (default: output.xlsx)')

    args = parser.parse_args()

    main(args)