import pysam
import pandas as pd
from collections import Counter
import numpy as np
from joblib import delayed,Parallel
import argparse
def read_is_split(read):
    if read.is_reverse:
        return read.cigartuples[0][0] == 4
    else:
        return read.cigartuples[-1][0] == 4
def get_split_site(read):
    if read.is_reverse:
        return read.reference_start
    else:
        return read.reference_end
def get_split_counter_and_strand(bamfile, row):
    split_counter = {True:Counter(), False: Counter()}
    supp_split_counter = Counter()
    strand_counter = Counter()
    bam = pysam.AlignmentFile(bamfile, 'rb')
    for read in bam.fetch(row['chrom'], row['start'], row['end']):
        if read_is_split(read) and read.mapping_quality > 0 and not read.is_supplementary:
            split_counter[read.is_reverse].update({get_split_site(read): 1})
        elif read.is_supplementary:
            supp_split_counter.update({get_split_site(read): 1})
        if not read.is_supplementary:
            strand_counter.update({read.is_reverse: 1})
    return split_counter, supp_split_counter,strand_counter
def get_softclip_sequence(read):
    if read.is_reverse:
        return read.query_sequence[:read.cigartuples[0][1]][::-1]
    else:
        return read.query_sequence[-read.cigartuples[-1][1]:]
def softclip_sequence_list(bamfile, row, loc, strand):
    sequence_list = []
    bam = pysam.AlignmentFile(bamfile, 'rb')
    for read in bam.fetch(row['chrom'], max(loc-10,0), loc+10):
        if read_is_split(read) and read.mapping_quality > 0 and not read.is_supplementary and strand == read.is_reverse:
            sequence_list.append(get_softclip_sequence(read))
    bam.close()
    return sequence_list
def build_consensus(l):
    max_len = max([len(s) for s in l])
    c = {i: Counter() for i in range(max_len)}
    for s in l:
        d = {i: {c:1} for i,c in enumerate(s)}
        for i, v in d.items():
            c[i].update(v)
    return ''.join([i_d.most_common()[0][0] for i, i_d in c.items()])
def make_insertion_stats(bamfile, row, idx):
    stats_dict = {}
    split_counter, supp_split_counter,strand_counter = get_split_counter_and_strand(bamfile, row)
    stats_dict['reads_pos'] = strand_counter[True]
    stats_dict['reads_neg'] = strand_counter[False]
    stats_dict['fraction_pos'] = strand_counter[True]/sum(strand_counter.values())
    stats_dict['split_reads_pos'] = sum(split_counter[True].values())
    stats_dict['split_reads_neg'] = sum(split_counter[False].values())
    if stats_dict['split_reads_pos'] > 0:
        stats_dict['most_common_pos_split'] = split_counter[True].most_common()[0][0]
        stats_dict['most_common_pos_split_n_reads'] = split_counter[True].most_common()[0][1]
        stats_dict['most_common_pos_split_fraction'] = stats_dict['most_common_pos_split_n_reads']/stats_dict['split_reads_pos']
        stats_dict['split_consensus_seq_pos'] = build_consensus(softclip_sequence_list(bamfile, row,stats_dict['most_common_pos_split'], True))
    else:
        stats_dict['most_common_pos_split'] = np.nan
        stats_dict['most_common_pos_split_n_reads'] = np.nan
        stats_dict['most_common_pos_split_fraction'] = np.nan
        stats_dict['split_consensus_seq_pos'] = ''
    if stats_dict['split_reads_neg'] > 0:
        stats_dict['most_common_neg_split'] = split_counter[False].most_common()[0][0]
        stats_dict['most_common_neg_split_n_reads'] = split_counter[False].most_common()[0][1]
        stats_dict['most_common_neg_split_fraction'] = stats_dict['most_common_neg_split_n_reads']/stats_dict['split_reads_neg']
        stats_dict['split_consensus_seq_neg'] = build_consensus(softclip_sequence_list(bamfile, row,stats_dict['most_common_neg_split'], False))
    else:
        stats_dict['most_common_neg_split'] = np.nan
        stats_dict['most_common_neg_split_n_reads'] = np.nan
        stats_dict['most_common_neg_split_fraction'] = np.nan
        stats_dict['split_consensus_seq_neg'] = ''
    return idx, stats_dict
def determine_strand(row):
    if row['fraction_pos'] > 0.9:
        return 'Positive'
    elif row['fraction_pos'] < 0.1:
        return 'Negative'
    else:
        return 'Ambiguous'
def make_pwm_df(df):
    l = df.values
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
    parser.add_argument('-bam','--bam_input',metavar='input', type=str, help='Input discordant .bam file')
    parser.add_argument('-bed','--bed_input', metavar='input', type=str, help='Input .bed file')
    parser.add_argument('-o','--output', metavar='output', type=str, help='Output bed file')
    parser.add_argument('-t','--threads', metavar='threads', type=int, help='threads')
    args = parser.parse_args()

    bamfile = args.bam_input
    bed_file = args.bed_input
    bed_output = args.output
    threads = args.threads

    bed_df = pd.read_csv(bed_file, sep='\t', header=None)
    bed_df.columns = ['chrom', 'start', 'end', 'n', 'mean_mapq','n_tag_sites']

    res = Parallel(n_jobs=threads, verbose = 3, backend='loky')(delayed(make_insertion_stats)(bamfile, row, idx) for idx, row in bed_df.iterrows())
    res_dict = {idx: stats_dict for idx, stats_dict in res}
    stats_df = pd.DataFrame.from_dict(res_dict, orient='index')
    stats_df = stats_df.join(bed_df)
    stats_df['strand'] = stats_df.apply(lambda row: determine_strand(row), axis=1)
    stats_df['total_reads'] = stats_df['reads_pos']+stats_df['reads_neg']

    split_read_sequence_count = {}
    for idx,row in stats_df.iterrows():
        split_read_sequence_count[idx] = {}
        split_read_sequence_count[idx]['seq_neg_len'] = len(row['split_consensus_seq_neg'])
        split_read_sequence_count[idx]['seq_pos_len'] = len(row['split_consensus_seq_pos'])
        split_read_sequence_count[idx]['seq_neg_T'] = row['split_consensus_seq_neg'][:10].count('T')
        split_read_sequence_count[idx]['seq_pos_A'] = row['split_consensus_seq_pos'][:10].count('A')
    split_read_stats = pd.DataFrame.from_dict(split_read_sequence_count, orient='index')
    stats_df = stats_df.join(split_read_stats)

    stats_df.to_csv(bed_output)

if __name__ == '__main__':
    main()
