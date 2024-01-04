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
    for read in bam_discord.fetch(row['chrom'],row['start'],row['end']):
        if not read_is_softclip(read) and  not read.is_supplementary and read.mapping_quality > 0:
            tag_site_counter.update({get_tag_site(read):1})
    return tag_site_counter

def filter_by_tag_sites(bed_in, bamfile, bed_out):
    bed_df = pd.read_csv(bed_in, sep='\t', header=None)
    bed_df.columns = ['chrom', 'start', 'end', 'n', 'mean_mapq']
    bed_df = bed_df[bed_df['mean_mapq'] > 0]
 
    bam_discord = pysam.AlignmentFile(bamfile, 'rb')
    
    tag_site_n_dict = {}
    for idx, row in bed_df.iterrows():
        tag_site_counter = get_tag_site_counter(row,bam_discord)
        tag_site_set = set([k for k,v in tag_site_counter.items() if v > 10])
        tag_site_n_dict[idx] = len(tag_site_set)
    tag_site_n_series = pd.Series(tag_site_n_dict)
    tag_site_n_series.name = 'n_tag_sites'
    bed_df = bed_df.join(tag_site_n_series)

    bed_df_filtered = bed_df[bed_df['n_tag_sites'] > 0]

    with open(bed_out, 'w') as f:
        for idx,row in bed_df_filtered.iterrows():
            f.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(row['chrom'], row['start'], row['end'], row['n'], row['mean_mapq'], row['n_tag_sites']))

def main():
    parser = argparse.ArgumentParser(description='Filter bed file', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i','--input',metavar='input', type=str, help='Input .bed file')
    parser.add_argument('-b','--input_discon', metavar='output', type=str, help='Input discordant .bam file')
    parser.add_argument('-o','--output', metavar='output', type=str, help='Output bed file')
    args = parser.parse_args()

    filter_by_tag_sites(args.input, args.input_discon, args.output)


if __name__ == '__main__':
    main()   
