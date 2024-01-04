import pysam
from collections import Counter
import pandas as pd
import argparse
from joblib import delayed,Parallel

def get_tag_site(read):
    if read.is_reverse:
        return read.reference_end
    else:
        return read.reference_start
def count_tag_sites_and_read_depth(bamfile, cell_BC_set, idx, row, tag):
    cell_counter_dict = {cell: Counter() for cell in cell_BC_set}
    bam = pysam.AlignmentFile(bamfile, 'rb')
    for read in bam.fetch(row['chrom'], row['start'], row['end']):
        cell = read.get_tag(tag)
        if read.mapping_quality > 0 and not read.is_supplementary:
            cell_counter_dict[cell].update({get_tag_site(read): 1})
    cell_counter_dict = {cell: c for cell,c in cell_counter_dict.items() if len(c) > 0}
    return idx, cell_counter_dict

def make_stats_per_donor(samplesheet, stats_df, read_depth_df, tag_site_df, prefix):
    for donor, donor_df in samplesheet.groupby('Donor'):
        donor_BCs = donor_df.index
        read_depth_df_subset = read_depth_df.reindex(donor_BCs, axis=1).fillna(0)
        tag_site_df_subset = tag_site_df.reindex(donor_BCs, axis=1).fillna(0)
        total_read_depth = read_depth_df_subset.sum(axis=1)
        any_detection = (read_depth_df_subset > 0).sum(axis=1)
        over_10_detection = (read_depth_df_subset > 10).sum(axis=1)
        over_1_detection = (read_depth_df_subset > 1).sum(axis=1)
        mean_tag_sites = tag_site_df_subset.mean(axis=1)

        total_read_depth.name = 'total_read_depth'
        any_detection.name = 'any_detection'
        over_10_detection.name = 'over_10_detection'
        over_1_detection.name = 'over_1_detection'
        mean_tag_sites.name = 'mean_tag_sites'

        additional_df = pd.DataFrame([total_read_depth, any_detection, over_10_detection, over_1_detection, mean_tag_sites]).T
        
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
    parser.add_argument('-bam-polyA', '--bam_input_polyA',  metavar='polyA_input', type='str', help='Remapped .bam file')
    parser.add_argument('-stats','--stats', metavar='stats', type=str, help='Input stats file')
    parser.add_argument('-sample','--samplesheet', metavar='samplesheet', type=str, help='Samplesheet file')
    parser.add_argument('-p','--prefix', metavar='prefix', type=str, help='prefix')
    parser.add_argument('-ct','--tag', metavar='tag', type=str, help='cell tag')
    parser.add_argument('-t','--threads', metavar='threads', type=int, help='threads')
    args = parser.parse_args()

    bamfile = args.bam_input
    bamfile_polyA = args.bam_input_polyA
    stats_file = args.stats
    samplesheet_file = args.samplesheet
    prefix = args.prefix
    tag = args.tag
    threads = args.threads
    samplesheet = pd.read_csv(samplesheet_file, index_col=0)
    cell_BC_set = set(samplesheet.index)

    stats_df = pd.read_csv(stats_file, index_col=0)

    res = Parallel(n_jobs=threads, verbose = 3, backend='loky')(delayed(count_tag_sites_and_read_depth)(bamfile, cell_BC_set, idx, row, tag) for idx, row in stats_df.iterrows())
    cell_counter_dict_total = {idx: cell_counter_dict for idx, cell_counter_dict in res}
    read_depth_df, tag_site_df = make_tables(cell_counter_dict_total,cell_BC_set)
    make_stats_per_donor(samplesheet, stats_df, read_depth_df, tag_site_df, prefix)

if __name__ == '__main__':
    main()
