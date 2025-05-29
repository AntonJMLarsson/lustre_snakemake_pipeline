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
    
def count_tag_sites_and_read_depth(bamfile, bamfile_polyA, cell_BC_set, idx, row, tag):
    tag_sites_per_cell_counter = {}
    
    bam = pysam.AlignmentFile(bamfile, 'rb')
    for read in bam.fetch(row['chrom'], row['start'], row['end']):
        cell = read.get_tag(tag)
        if cell not in tag_sites_per_cell_counter and cell in cell_BC_set:
            tag_sites_per_cell_counter[cell] = Counter()
        if read.mapping_quality > 0 and not read.is_supplementary and not read_is_softclip(read):
            tag_sites_per_cell_counter[cell].update({get_tag_site(read): 1})
    bam.close()

    bam = pysam.AlignmentFile(bamfile_polyA, 'rb')
    for read in bam.fetch(row['chrom'], row['start'], row['end']):
        cell = read.get_tag(tag)
        if cell not in tag_sites_per_cell_counter and cell in cell_BC_set:
            tag_sites_per_cell_counter[cell] = Counter()
        if read.mapping_quality > 0 and not read.is_supplementary and not read_is_softclip(read):
            tag_sites_per_cell_counter[cell].update({get_tag_site(read): 1})
    bam.close()

    res = []
    for cell, cell_tag_site_counter in tag_sites_per_cell_counter.items():
        for tag_site, count in cell_tag_site_counter.items():
            res.append([cell, idx, tag_site, count])
    return res

def main():
    parser = argparse.ArgumentParser(description='Count cell detection', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-bam','--bam_input',metavar='input', type=str, help='Input discordant .bam file')
    parser.add_argument('-bam-polyA', '--bam_input_polyA',  metavar='polyA_input', type=str, help='Remapped .bam file')
    parser.add_argument('-stats','--stats', metavar='stats', type=str, help='Input stats file')
    parser.add_argument('-sample','--samplesheet', metavar='samplesheet', type=str, help='Samplesheet file')
    parser.add_argument('-d','--donor', metavar='donor', type=str, help='Donor')
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
    cell_BC_set = set(samplesheet[samplesheet['Donor'] == args.donor].index)

    stats_df = pd.read_csv(stats_file, index_col=0)

    res = Parallel(n_jobs=threads, verbose = 3, backend='loky')(delayed(count_tag_sites_and_read_depth)(bamfile, bamfile_polyA, cell_BC_set, idx, row, tag) for idx, row in stats_df.iterrows())
    cell_counter_dict_total = {idx: cell_counter_dict for idx, cell_counter_dict in res}
    read_depth_df, tag_site_df = make_tables(cell_counter_dict_total,cell_BC_set)
    make_stats_per_donor(samplesheet, stats_df, read_depth_df, tag_site_df, prefix)

if __name__ == '__main__':
    main()
