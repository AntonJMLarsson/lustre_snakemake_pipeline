import pandas as pd
import argparse
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
def filter_stats_file(stats_file,prefix):
    df = pd.read_csv(stats_file, index_col=0)
    df['pass'] = df.apply(lambda row: is_passing(row), axis=1)
    df_filtered = df[df['pass']]
    df_filtered = df_filtered[df_filtered['total_read_depth']/df['total_reads'] > 0.2]
    df_filtered.to_csv('{}_stats_filtered.csv'.format(prefix))

    with open('{}_stats_file.bed'.format(prefix), 'w') as f:
        for idx, row in df[df['pass']].iterrows():
            f.write('{}\t{}\t{}\t{}\n'.format(row['chrom'], row['start'], row['end'], '{}_{}'.format(prefix,idx)))
def main():
    parser = argparse.ArgumentParser(description='Count cell detection', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-stats','--stats', metavar='stats', type=str, help='Input stats file')
    parser.add_argument('-p','--prefix', metavar='prefix', type=str, help='prefix')

    args = parser.parse_args()

    stats_file = args.stats
    prefix = args.prefix

    filter_stats_file(stats_file, prefix)
    
if __name__ == '__main__':
    main()
