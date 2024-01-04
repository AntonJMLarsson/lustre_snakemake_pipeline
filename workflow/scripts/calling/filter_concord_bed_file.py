import pysam
import pandas as pd
import argparse
def read_is_not_supplementary(read):
    return not read.is_supplementary
def read_is_not_supplementary_concord(read):
    return not read.is_supplementary and read.is_proper_pair

def filter_file(bed_input, bamfile_concord, bamfile_discord, bed_output):
    df_bed = pd.read_csv(bed_input, sep='\t', header=None)
    df_bed.columns = ['chrom', 'start', 'end', 'mean_mapq']
    bam_concord = pysam.AlignmentFile(bamfile_concord, 'rb')
    bam_discord = pysam.AlignmentFile(bamfile_discord, 'rb')
    
    count_dict = {}
    for i, row in df_bed.iterrows():
        concord_count = bam_concord.count(row['chrom'], row['start'], row['end'], read_callback=read_is_not_supplementary_concord)
        discord_count = bam_discord.count(row['chrom'], row['start'], row['end'], read_callback=read_is_not_supplementary)
        count_dict[i] = (concord_count, discord_count, concord_count+discord_count)
    df = pd.DataFrame(count_dict, index=['concord_count', 'discord_count', 'total_count']).T
    df['fraction_concord'] = df['concord_count']/df['total_count']
    with open(bed_output, 'w') as f:
        for idx in df[df['fraction_concord'] > 0.01].index:
            row = df_bed.loc[idx]
            f.write('{}\t{}\t{}\t{}\n'.format(row['chrom'], row['start'], row['end'], row['mean_mapq']))
    return None
def main():
    parser = argparse.ArgumentParser(description='Filter bed file', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i','--input',metavar='input', type=str, help='Input .bed file')
    parser.add_argument('-id','--input_discon', metavar='output', type=str, help='Input discordant .bam file')
    parser.add_argument('-ic','--input_concord', metavar='output', type=str, help='Input concordant .bam file')
    parser.add_argument('-o1','--output', metavar='output', type=str, help='Output bed file')
    args = parser.parse_args()
    

    filter_file(args.input, args.input_concord, args.input_discon, args.output) 

if __name__ == '__main__':
    main()
