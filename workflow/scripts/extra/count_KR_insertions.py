import pysam
from collections import Counter
import pandas as pd
import argparse
from joblib import delayed,Parallel

def count_tag_sites_and_read_depth(bamfile, cell_BC_set, idx, row, tag):
    cell_counter_dict = Counter()
    bam = pysam.AlignmentFile(bamfile, 'rb')
    if row['strand'] == '+':
        start = row['start']
        end = row['end']+5000
    else:
        start = row['start']-5000
        end = row['end']
    try:
        for read in bam.fetch(row['chrom'], start, end):
            cell = read.get_tag(tag)
            if read.mapping_quality > 0 and not read.is_supplementary:
                cell_counter_dict.update({cell: 1})
    except:
        return idx, None
    return idx, cell_counter_dict
def make_tables(cell_counter_dict_total):
    read_depth_dict = {}
    for idx, cell_counter_dict in cell_counter_dict_total.items():
        read_depth_dict[idx] = {}
        for cell, cell_count in cell_counter_dict.items():
            read_depth_dict[idx][cell] = cell_count
    read_depth_df = pd.DataFrame.from_dict(read_depth_dict,orient='index').fillna(0)
    return read_depth_df

def make_stats_per_donor(samplesheet, stats_df, read_depth_df, prefix):
    for donor, donor_df in samplesheet.groupby('Donor'):
        donor_BCs = donor_df.index
        read_depth_df_subset = read_depth_df.reindex(donor_BCs, axis=1).fillna(0)
        total_read_depth = read_depth_df_subset.sum(axis=1)
        any_detection = (read_depth_df_subset > 0).sum(axis=1)
        over_10_detection = (read_depth_df_subset > 10).sum(axis=1)
        

        total_read_depth.name = 'total_read_depth'
        any_detection.name = 'any_detection'
        over_10_detection.name = 'over_10_detection'

        additional_df = pd.DataFrame([total_read_depth, any_detection, over_10_detection]).T

        stats_df_donor = stats_df.join(additional_df)


        read_depth_df_subset.to_csv('{}_{}_read_depth.csv'.format(prefix, donor))
        stats_df_donor.to_csv('{}_{}_stats.csv'.format(prefix, donor))
    return True

def main():
    parser = argparse.ArgumentParser(description='Count KR cell detection', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-bam','--bam_input',metavar='input', type=str, help='Input discordant .bam file')
    parser.add_argument('-ins','--insertions', metavar='stats', type=str, help='Input insertion bed file')
    parser.add_argument('-sample','--samplesheet', metavar='samplesheet', type=str, help='Samplesheet file')
    parser.add_argument('-p','--prefix', metavar='prefix', type=str, help='prefix')
    parser.add_argument('-ct','--tag', metavar='tag', type=str, help='cell tag')
    parser.add_argument('-t','--threads', metavar='threads', type=int, help='threads')
    args = parser.parse_args()

    bamfile = args.bam_input
    bed_file = args.insertions
    samplesheet_file = args.samplesheet
    prefix = args.prefix
    tag = args.tag
    threads = args.threads
    samplesheet = pd.read_csv(samplesheet_file, index_col=0)
    cell_BC_set = set(samplesheet.index)

    bed_df = pd.read_csv(bed_file, sep='\t', header=None)
    bed_df.columns = ['chrom', 'start', 'end', 'name', 'bla', 'strand']

    res = Parallel(n_jobs=threads, verbose = 3, backend='loky')(delayed(count_tag_sites_and_read_depth)(bamfile, cell_BC_set, idx, row, tag) for idx, row in bed_df.iterrows())

    cell_counter_dict_total = {idx: cell_counter_dict for idx, cell_counter_dict in res if cell_counter_dict is not None}
    read_depth_df = make_tables(cell_counter_dict_total)

    make_stats_per_donor(samplesheet, bed_df, read_depth_df, prefix)

if __name__ == '__main__':
    main()