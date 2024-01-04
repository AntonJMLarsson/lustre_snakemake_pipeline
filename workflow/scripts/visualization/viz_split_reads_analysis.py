import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
def make_pwm_df(df):
    l = df.values
    l = [s for s in l if type(s) == str]
    max_len = max([len(s) for s in l])
    c = {i: Counter() for i in range(max_len)}
    for s in l:
        d = {i: {c:1} for i,c in enumerate(s)}
        for i, v in d.items():
            c[i].update(v)
    pwm = {}
    for i,v in c.items():
        pwm[i] = {}
        s = sum(v.values())
        for char, count in v.items():
            pwm[i][char] = count/s
    pwm_df = pd.DataFrame(pwm)
    pwm_df = pwm_df.fillna(0)
    return pwm_df
def main():
    parser = argparse.ArgumentParser(description='Analyze split reads', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i','--input', metavar='input', type=str, help='Input .csv file')
    parser.add_argument('-p','--prefix', metavar='output', type=str, help='Output prefix')
    parser.add_argument('-q','--mapq_filter', metavar='mapq', type=int, default=25, help='Only plot insertions with mean mapq > q')
    args = parser.parse_args()

    input_file = args.input
    prefix = args.prefix
    mapq_filter = args.mapq_filter

    stats_df = pd.read_csv(input_file, index_col = 0)
    stats_df_over_mapq = stats_df[stats_df['mean_mapq'] > mapq_filter]

    pwm_df_dict = {}
    for strand, stats_df_strand in stats_df_over_mapq.groupby('strand'):
        pwm_df_dict[strand] = {}
        pwm_df_dict[strand]['neg'] = make_pwm_df(stats_df_strand['split_consensus_seq_neg']).T
        pwm_df_dict[strand]['pos'] = make_pwm_df(stats_df_strand['split_consensus_seq_pos']).T
    color_dict = {'A': 'Green', 'C': 'Blue', 'T': 'Red', 'G': 'Yellow'}
    group_to_idx = {'Positive': 0, 'Negative': 1, 'Ambiguous': 2}
    strand_to_idx = {'neg': 1, 'pos':0}
    fig, axes = plt.subplots(3,2, figsize=(10,10), sharey=True)
    for group, new_dict in pwm_df_dict.items():
        for strand, new_df in new_dict.items():
            new_df.plot(ax=axes[group_to_idx[group]][strand_to_idx[strand]],color=[color_dict.get(x, '#333333') for x in new_df.columns])
    axes[0][0].set_title('Split positive direction')
    axes[0][1].set_title('Split negative direction')
    axes[0][0].set_ylabel('Positive insertions \n fraction of bases')
    axes[1][0].set_ylabel('Negative insertions \n fraction of bases')
    axes[2][0].set_ylabel('Ambiguous insertions \n fraction of bases')
    axes[2][0].set_xlabel('Position after split (bp)')
    axes[2][1].set_xlabel('Position after split (bp)')
    plt.tight_layout()
    plt.savefig('{}_split_reads_sequences.pdf'.format(prefix))
    plt.clf()

    fig, ax = plt.subplots(1,2)

    sns.histplot(hue='strand', x='seq_pos_A', data=stats_df_over_mapq[stats_df_over_mapq['seq_pos_len'] > 9],ax=ax[0])
    sns.histplot(hue='strand', x='seq_neg_T', data=stats_df_over_mapq[stats_df_over_mapq['seq_neg_len'] > 9],ax=ax[1])
    ax[0].set_xlabel('Number of As within first 10 bases')
    ax[1].set_xlabel('Number of Ts within first 10 bases')
    ax[0].set_title('Split reads positive direction')
    ax[1].set_title('Split reads negative direction')
    plt.tight_layout()
    plt.savefig('{}_first_10_split_bases.pdf'.format(prefix))
    plt.clf()

    plt.hist(stats_df_over_mapq['fraction_pos'], bins=30)
    plt.xlabel('Fraction of reads \n supporting positive insertion')
    plt.ylabel('Count (Insertions)')
    plt.tight_layout()
    plt.savefig('{}_strand_assignment.pdf'.format(prefix))
    plt.clf()

    fig, ax = plt.subplots(1,2, sharey=True)

    sns.stripplot(y='most_common_pos_split_fraction', x='strand',data=stats_df_over_mapq,ax=ax[0])
    sns.stripplot(y='most_common_neg_split_fraction', x='strand', data=stats_df_over_mapq,ax=ax[1])
    ax[0].set_ylabel('Fraction support for split position')
    ax[1].set_ylabel('Fraction support for split position')
    ax[0].set_title('Split reads positive direction')
    ax[1].set_title('Split reads negative direction')
    plt.tight_layout()
    plt.savefig('{}_split_support.pdf'.format(prefix))
    plt.clf()

if __name__ == "__main__":
    main()

