import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse

def any_split_reads(row):
    if row['fraction_pos'] > 0.5:
        if row['split_reads_pos'] > 0 or row['split_reads_pos_remapped'] > 0:
            return True
        else:
            return False
    else:
        if row['split_reads_neg'] > 0 or row['split_reads_neg_remapped'] > 0:
            return True
        else:
            return False

def is_passing(row):
    if row['mean_mapq'] > 25 and row['over_10_detection'] > 0:
        if row['strand'] == 'Positive' or row['fraction_pos'] > 0.5:
            if row['seq_pos_len'] > 9:
                if row['seq_pos_A'] > 7:
                    return True
                else:
                    return False
            elif row['seq_pos_len'] > -1:
                if row['seq_pos_len'] == row['seq_pos_A']:
                    return True
                else:
                    return False
            else:
                return False
        elif row['strand'] == 'Negative' or row['fraction_pos'] < 0.5:
            if row['seq_neg_len'] > 9:
                if row['seq_neg_T'] > 7:
                    return True
                else:
                    return False
            elif row['seq_neg_len'] > -1:
                if row['seq_neg_len'] == row['seq_neg_T']:
                    return True
                else:
                    return False
            else:
                return False
        else:
            False
    else:
        return False
    return False

def apply_filtering_steps(regular_stats, remapped_stats, read_matrix, read_matrix_remapped, prefix, name):
    regular_df = pd.read_csv(regular_stats, index_col=0)
    regular_df.index = regular_df.apply(lambda row: '{}-{}'.format(name, row.name), axis=1)
    remapped_df = pd.read_csv(remapped_stats, index_col=0)
    remapped_df.index = remapped_df.apply(lambda row: '{}-{}'.format(name, row.name), axis=1)
    read_mat = pd.read_csv(read_matrix, index_col=0)
    read_mat.index = read_mat.apply(lambda row: '{}-{}'.format(name, row.name), axis=1)
    read_mat.columns = read_mat.apply(lambda col: '{}-{}'.format(col.name, name))
    read_mat_remapped = pd.read_csv(read_matrix_remapped, index_col=0)
    read_mat_remapped.index = read_mat_remapped.apply(lambda row: '{}-{}'.format(name, row.name), axis=1)
    read_mat_remapped.columns = read_mat_remapped.apply(lambda col: '{}-{}'.format(col.name, name))

    full_df = regular_df.join(remapped_df, rsuffix='_remapped')
    full_df = full_df[full_df['mean_mapq'] > 25]

    full_df['pass'] = full_df.apply(lambda row: is_passing(row), axis=1)
    full_df_filtered = full_df[full_df['pass']]

    full_df_discarded = full_df[~full_df['pass']]

    full_df_filtered_1_tag_site = full_df_filtered[(full_df_filtered['n_tag_sites'] == 1)][full_df_filtered[full_df_filtered['n_tag_sites'] == 1].apply(lambda row: any_split_reads(row), axis=1)]
    full_df_conservative = full_df_filtered[full_df_filtered['n_tag_sites'] > 1][full_df_filtered[full_df_filtered['n_tag_sites'] > 1].apply(lambda row: any_split_reads(row), axis=1)]
    full_df_conservative = full_df_conservative[full_df_conservative['total_read_depth']/full_df_conservative['total_reads'] > 0.2]
    full_df_both = full_df_conservative.append(full_df_filtered_1_tag_site)

    read_mat_both = read_mat + read_mat_remapped
    
    read_mat_both = read_mat_both.reindex(full_df_both.index)

    read_mat_both.to_csv('{}_read_matrix.csv'.format(prefix))

    full_df_both.to_csv('{}_insertion_table.csv'.format(prefix))

    full_df_discarded.to_csv('{}_discarded_insertions.csv'.format(prefix))

    return True



def main():
    parser = argparse.ArgumentParser(description='Prepare insertion table')
    parser.add_argument('--prefix', type=str, help='Prefix')
    parser.add_argument('--name', type=str, help='Experiment name')
    parser.add_argument('--regular', type=str, help='File of regular stats')
    parser.add_argument('--remapped', type=str, help='File of remapped stats')
    parser.add_argument('--read-mat', type=str, help='Read count matrix')
    parser.add_argument('--read-mat-remapped', type=str, help='Read count matrix remapped')
    

    args = parser.parse_args()

    prefix = args.prefix
    regular_stats = args.regular
    remapped_stats = args.remapped
    read_matrix = args.read_mat
    read_matrix_remapped = args.read_mat_remapped
    name = args.name

    apply_filtering_steps(regular_stats, remapped_stats, read_matrix, read_matrix_remapped, prefix, name)


if __name__ == '__main__':
    main()
