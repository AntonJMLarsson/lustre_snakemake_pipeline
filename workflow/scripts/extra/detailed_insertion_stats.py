import pysam
from collections import Counter
import pandas as pd
import argparse
from joblib import delayed,Parallel

def read_is_softclip(read):
    if read.is_reverse:
        return read.cigartuples[-1][0] == 4
    else:
        return read.cigartuples[0][0] == 4
    
def read_has_deletion(read):
    return 'D' in read.cigarstring

def get_tag_site(read):
    if read.is_reverse:
        return read.reference_end
    else:
        return read.reference_start

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
