import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import argparse

def make_plot(file, group, prefix):

    df = pd.read_csv(file, index_col=0)

    primer_columns = ['primer_distance_{}'.format(i) for i in range(24)]
    continue_columns = ['continue_distance_{}'.format(i) for i in range(9)]

    pass_primer = df[['primer_distance_{}'.format(i) for i in range(7)]].sum(axis=1)
    pass_primer.name = 'pass_primer'

    pass_continue = df[['continue_distance_{}'.format(i) for i in range(2)]].sum(axis=1)
    pass_continue.name = 'pass_continue'

    df = df.join(pass_primer)
    df = df.join(pass_continue)

    df['fraction_pass_primer'] = df['pass_primer']/df['total_n']

    df['fraction_pass_continue'] = df['pass_continue']/df['total_n']


    if group is None:
        group = 'Sample'
        df['Sample'] = 'ALL'

    primer_df = df[primer_columns]
    continue_df = df[continue_columns]
    continue_df_norm = continue_df.div(continue_df.sum(axis=1), axis=0).join(df[group])
    primer_df_norm = primer_df.div(primer_df.sum(axis=1), axis=0)
    n_groups = primer_df_norm.join(df[group]).groupby(group).count().shape[0]
    fig, axes = plt.subplots(n_groups,3, figsize=(16,2*n_groups), sharey=True, squeeze=False)
    fraction_pass_primer = primer_df_norm[['primer_distance_{}'.format(i) for i in range(7)]].sum(axis=1)
    i = 0
    for s, df_subset in primer_df_norm.join(df[group]).groupby(group):
        g1 = sns.boxplot(data=df_subset.drop(group, axis=1), ax=axes[i][0], showfliers=False, **{'boxprops':{'facecolor':'none'}})
        g1.set(xticklabels=[])
        g1.set(title=s)
        g1.set(ylim=(0,1))
        continue_subset = continue_df_norm[continue_df_norm[group] == s]
        g2 = sns.boxplot(data=continue_subset.drop(group, axis=1), ax=axes[i][1], showfliers=False, **{'boxprops':{'facecolor':'none'}})
        g2.set(xticklabels=[])
        g2.set(title=s)
        fraction_df = df[df[group] == s][['fraction_pass_primer', 'fraction_pass_continue', 'fraction_pass']]
        g3 = sns.boxplot(data=fraction_df, ax=axes[i][2])
        g3.set(xticklabels=[])
        g3.set(title=s)
        i += 1
        if i == n_groups:
            g1.set(xticklabels=range(24))
            g2.set(xticklabels=range(9))
            g1.set(xlabel='Distance from primer')
            g2.set(xlabel='Distance from continuation')
            g3.set(xticklabels=['Pass primer', 'Pass continue', 'Pass total'])
    plt.tight_layout()
    plt.savefig('{}_passing_reads_plot.pdf'.format(prefix))
    plt.show()

def main():
    parser = argparse.ArgumentParser(description='Analyze split reads', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i','--input', metavar='input', type=str, help='Input .csv file')
    parser.add_argument('-g', '--group', type=str, default=None, help='group')
    parser.add_argument('-p','--prefix', metavar='output', type=str, help='Output prefix')

    args = parser.parse_args()

    make_plot(args.input, args.group, args.prefix)

if __name__ == '__main__':
    main()