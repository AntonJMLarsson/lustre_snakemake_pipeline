import pysam
from collections import Counter
import pandas as pd
import numpy as np
import argparse
from joblib import delayed,Parallel
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
    for read in bam.fetch(row['chrom'], max(row['start'],0), row['end']+500):
        if read_is_split(read) and read.mapping_quality > 0 and not read.is_supplementary:
            split_counter[read.is_reverse].update({get_split_site(read): 1})
        elif read.is_supplementary:
            supp_split_counter.update({get_split_site(read): 1})
        if not read.is_supplementary:
            strand_counter.update({read.is_reverse: 1})
    bam.close()
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
def make_insertion_stats_remapped(bamfile, row, idx):
    stats_dict = {}
    split_counter, supp_split_counter,strand_counter = get_split_counter_and_strand(bamfile, row)
    stats_dict['reads_pos'] = strand_counter[True]
    stats_dict['reads_neg'] = strand_counter[False]
    if sum(strand_counter.values()) == 0:
        stats_dict['fraction_pos'] = np.nan
    else:
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

def get_tag_site(read):
    if read.is_reverse:
        return read.reference_end
    else:
        return read.reference_start
def count_tag_sites_remapped_reads(bamfile, cell_BC_set, idx, row, tag):
    cell_counter_dict = {cell: Counter() for cell in cell_BC_set}
    bam = pysam.AlignmentFile(bamfile, 'rb')
    for read in bam.fetch(row['chrom'], max(row['start']-500,0), row['end']+500):
        cell = read.get_tag(tag)
        if read.mapping_quality > 0 and not read.is_supplementary and read_is_split(read):
            cell_counter_dict[cell].update({get_tag_site(read): 1})
    cell_counter_dict = {cell: c for cell,c in cell_counter_dict.items() if len(c) > 0}
    bam.close()
    return idx, cell_counter_dict
def make_stats_per_donor(samplesheet, stats_df, read_depth_df, tag_site_df, prefix):
    for donor, donor_df in samplesheet.groupby('Donor'):
        donor_BCs = donor_df.index
        read_depth_df_subset = read_depth_df.reindex(donor_BCs, axis=1).fillna(0)
        tag_site_df_subset = tag_site_df.reindex(donor_BCs, axis=1).fillna(0)
        total_read_depth = read_depth_df_subset.sum(axis=1)
        any_detection = (read_depth_df_subset > 0).sum(axis=1)
        over_10_detection = (read_depth_df_subset > 10).sum(axis=1)
        mean_tag_sites = tag_site_df_subset.mean(axis=1)

        total_read_depth.name = 'total_read_depth'
        any_detection.name = 'any_detection'
        over_10_detection.name = 'over_10_detection'
        mean_tag_sites.name = 'mean_tag_sites'

        additional_df = pd.DataFrame([total_read_depth, any_detection, over_10_detection, mean_tag_sites]).T
        
        stats_df_donor = stats_df.join(additional_df)
        
        
        read_depth_df_subset.to_csv('{}_{}_read_depth.csv'.format(prefix, donor))
        tag_site_df_subset.to_csv('{}_{}_tag_sites.csv'.format(prefix, donor))
        stats_df_donor.to_csv('{}_{}_stats.csv'.format(prefix, donor))
    return True


def make_tables(cell_counter_dict_total, cell_BC_set):
    read_depth_dict = {}
    tag_site_dict = {}
    for idx, cell_counter_dict in cell_counter_dict_total.items():
        read_depth_dict[idx] = {}
        tag_site_dict[idx] = {}
        for cell in cell_BC_set:
            if cell in cell_counter_dict:
                cell_counter = cell_counter_dict[cell]
            else:
                cell_counter = Counter()
            read_depth_dict[idx][cell] = sum(cell_counter.values())
            tag_site_dict[idx][cell] = len(cell_counter.keys())
    read_depth_df = pd.DataFrame.from_dict(read_depth_dict,orient='index')
    tag_site_df = pd.DataFrame.from_dict(tag_site_dict,orient='index')
    return read_depth_df, tag_site_df

def main():
    parser = argparse.ArgumentParser(description='Count cell detection', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-bam','--bam_input',metavar='input', type=str, help='Input discordant .bam file')
    parser.add_argument('-stats','--stats', metavar='stats', type=str, help='Input stats file')
    parser.add_argument('-sample','--samplesheet', metavar='samplesheet', type=str, help='Samplesheet file')
    parser.add_argument('-p','--prefix', metavar='prefix', type=str, help='prefix')
    parser.add_argument('-ct','--tag', metavar='tag', type=str, help='cell tag')
    parser.add_argument('-t','--threads', metavar='threads', type=int, help='threads')
    args = parser.parse_args()

    bamfile = args.bam_input
    stats_file = args.stats
    samplesheet_file = args.samplesheet
    prefix = args.prefix
    tag = args.tag
    threads = args.threads
    samplesheet = pd.read_csv(samplesheet_file, index_col=0)
    cell_BC_set = set(samplesheet.index)

    stats_df_old = pd.read_csv(stats_file, index_col=0)

    res = Parallel(n_jobs=threads, verbose = 3, backend='loky')(delayed(make_insertion_stats_remapped)(bamfile, row, idx) for idx, row in stats_df_old.iterrows())
    res_dict = {idx: stats_dict for idx, stats_dict in res}
    stats_df = pd.DataFrame.from_dict(res_dict, orient='index')
    stats_df['strand'] = stats_df.apply(lambda row: determine_strand(row), axis=1)
    stats_df['total_reads'] = stats_df['reads_pos']+stats_df['reads_neg']

    stats_df_full = stats_df_old.join(stats_df, rsuffix = '_remapped')

    res = Parallel(n_jobs=threads, verbose = 3, backend='loky')(delayed(count_tag_sites_remapped_reads)(bamfile, cell_BC_set, idx, row, tag) for idx, row in stats_df_full.iterrows())
    cell_counter_dict_total = {idx: cell_counter_dict for idx, cell_counter_dict in res}
    read_depth_df, tag_site_df = make_tables(cell_counter_dict_total,cell_BC_set)
    make_stats_per_donor(samplesheet, stats_df, read_depth_df, tag_site_df, prefix)


if __name__ == '__main__':
    main()
