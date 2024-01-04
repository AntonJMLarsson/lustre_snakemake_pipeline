import pysam
import pandas as pd
from collections import Counter
import argparse
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
def get_tag_site_counter(row,bam_discord):
    tag_site_counter = Counter()
    for read in bam_discord.fetch(row['chrom'],max(row['start']-500,0),row['end']+500):
        if not read_is_softclip(read) and not read.is_supplementary and read.mapping_quality > 0:
            tag_site_counter.update({get_tag_site(read):1})
    return tag_site_counter

def count_tag_sites(stats_in, bamfile_regular, bamfile_remapped, stats_out):
    bam_regular = pysam.AlignmentFile(bamfile_regular, 'rb')
    bam_remapped = pysam.AlignmentFile(bamfile_remapped, 'rb')

    df_stats = pd.read_csv(stats_in, index_col=0)

    tag_site_n_dict = {}
    for name, row in df_stats.iterrows():
        tag_site_counter_regular = get_tag_site_counter(row, bam_regular)
        tag_site_counter_remapped = get_tag_site_counter(row, bam_remapped)
        tag_site_counter_regular.update(tag_site_counter_remapped)

        tag_site_set = set([k for k,v in tag_site_counter_regular.items() if v > 4])

        tag_site_n_dict[name] = len(tag_site_set)
    tag_site_n_series = pd.Series(tag_site_n_dict)
    tag_site_n_series.name = 'n_tag_sites_w_remapped'

    df_stats = df_stats.join(tag_site_n_series)

    df_stats.to_csv(stats_out)

def main():
    parser = argparse.ArgumentParser(description='Count tag sites after remapping', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-bam','--bam_input',metavar='input', type=str, help='Input discordant .bam file')
    parser.add_argument('-bam-polyA', '--bam_input_polyA',  metavar='polyA_input', type=str, help='Remapped .bam file')
    parser.add_argument('-stats_in','--stats_in', metavar='stats', type=str, help='Input stats file')
    parser.add_argument('-stats_out','--stats_out', metavar='stats', type=str, help='Output stats file')

    args = parser.parse_args()

    stats_in = args.stats_in
    stats_out = args.stats_out
    bamfile_regular = args.bam_input
    bamfile_remapped = args.bam_input_polyA

    count_tag_sites(stats_in, bamfile_regular, bamfile_remapped, stats_out)
if __name__ == "__main__":
    main()


