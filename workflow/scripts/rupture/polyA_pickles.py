import argparse
import pickle
import glob
import statistics
from pandas import pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Process polyA pickle files and export per-cell summary table."
    )
    parser.add_argument(
        "insertion_csv",
        help="Path to insertion CSV file."
    )
    parser.add_argument(
        "samplesheet_csv",
        help="Path to samplesheet CSV file."
    )
    parser.add_argument(
        "output_csv",
        help="Path to output CSV file."
    )
    return parser.parse_args()


def main():
    args = parse_args()


    somatic_files = glob.glob('polyA_pkls/*.pkl')
    somatic_res = {}
    for filename in somatic_files:
        ins_name = '_'.join(filename.split('/')[1].split('_')[:-1])
        with open(filename, 'rb') as f:
            polyA_res = pickle.load(f)
            somatic_res[ins_name] = polyA_res


    insertion_df = pd.read_csv(args.insertion_csv, index_col=0)

    samplesheet_df = pd.read_csv(args.samplesheet_csv, index_col=0)

    if 'ins_name' not in insertion_df.columns:
        insertion_df['ins_name'] = insertion_df.apply(lambda row: '{}_{}_{}'.format(row['chrom'], row['start'], row['end']), axis=1)

    somatic_median_per_cell = {}
    read_count_per_cell = {}
    for ins_name, polyA_dict in somatic_res.items():
        somatic_median_per_cell[ins_name] = {}
        read_count_per_cell[ins_name] = {}
        for cell, counter in polyA_dict.items():
            somatic_median_per_cell[ins_name][cell] = statistics.median(counter)
            read_count_per_cell[ins_name][cell] = sum(counter.values())
        if len(somatic_median_per_cell[ins_name]) == 0:
            del somatic_median_per_cell[ins_name]
            del read_count_per_cell[ins_name]

    df_list = []
    for ins_name, median_dict in somatic_median_per_cell.items():
        if ins_name not in insertion_df['ins_name'].values:
            continue
        read_count_df = pd.DataFrame(read_count_per_cell[ins_name], index=['n_counts']).T
        df = pd.DataFrame(median_dict, index=['length']).T.join(samplesheet_df).join(read_count_df)
        df = df[df['n_counts'] > 2]
        if df.shape[0] == 0:
            continue
        df['ins_name'] = ins_name
        df_list.append(df)

    df_full = pd.concat(df_list)
    df_full.to_csv(args.output_csv)


if __name__ == '__main__':
    main()