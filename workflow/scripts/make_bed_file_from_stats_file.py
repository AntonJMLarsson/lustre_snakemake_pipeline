import pandas as pd
import os
import argparse

def make_bed_file_from_stats_df(file_list, bed_file):
    df_list = []
    for stats_file in file_list:
        df_stats = pd.read_csv(stats_file, index_col=0)
        name = os.path.splitext(stats_file)[0].split('/')[-1]
        df_stats['name'] = df_stats.apply(lambda row: '{}-{}'.format(name, row.name), axis=1)
        df_stats['score'] = 0
        df_stats['strand_for_bed'] = df_stats.apply(lambda row: '+' if row['fraction_pos'] > 0.5 else '-', axis=1)
        df_stats_bed = df_stats[['chrom', 'start', 'end', 'name', 'score', 'strand_for_bed']]
        df_list.append(df_stats_bed)
    df_full = pd.concat(df_list)
    df_full.to_csv(bed_file, sep='\t', header=None, index=False)

def main():
    parser = argparse.ArgumentParser(description='Make bed file', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i','--input', nargs='+', metavar='list', required=True, help='List of input files')
    parser.add_argument('-o','--output', metavar='bed', type=str, required=True, help='Output bed file')
    args = parser.parse_args()

    file_list = args.input
    bed_file = args.output

    make_bed_file_from_stats_df(file_list, bed_file)

if __name__ == "__main__":
    main()
