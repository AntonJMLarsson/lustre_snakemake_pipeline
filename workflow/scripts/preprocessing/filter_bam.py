import pysam
import argparse
from collections import Counter
import pandas as pd
L1_primer_seq = 'TGCACATGTACCCTAAAACTTAG'
def hamming_distance(chaine1, chaine2):
    return sum(c1 != c2 for c1, c2 in zip(chaine1, chaine2))

def filter_bam(bamfile, outbam, report, sample):
    bam_in = pysam.AlignmentFile(bamfile, 'rb')
    bam_out = pysam.AlignmentFile(outbam, 'wb', template=bam_in)
    prev_qname = ''
    read_list = []
    read2_stats = {}
    
    n_pass = 0
    n_fail = 0

    primer_distance_counter = Counter()
    continue_distance_counter = Counter()
    for read in bam_in.fetch(until_eof=True):
        if read.query_name != prev_qname:
            if len(read2_stats) != 0:
                primer_distance_counter.update({read2_stats['primer_distance']:1})
                continue_distance_counter.update({read2_stats['continue_distance']:1})
                if read2_stats['primer_distance'] < 7 and read2_stats['continue_distance'] < 2:
                    n_pass += 1
                    for quality_read in read_list:
                        bam_out.write(quality_read)
                else:
                    n_fail += 1
            else:
                n_fail += 1
            read_list = []
            read2_stats = {}
        if read.is_read2 and not read.is_supplementary and not read.is_unmapped:
            forward_seq = read.get_forward_sequence()
            read2_stats['primer_distance'] = hamming_distance(L1_primer_seq, forward_seq[:23])
            read2_stats['continue_distance'] = hamming_distance('AGTATAAT', forward_seq[23:31])
        read_list.append(read)
        prev_qname = read.query_name
    if len(read2_stats) != 0:
        if read2_stats['primer_distance'] < 7 and read2_stats['continue_distance'] < 2:
            n_pass += 1
            for quality_read in read_list:
                bam_out.write(quality_read)
        else:
            n_fail += 1

    full_df = pd.DataFrame([[n_pass, n_fail]], index=[sample], columns=['n_pass', 'n_fail'])
    primer_dist_df = pd.DataFrame(primer_distance_counter, index=[sample]).sort_index(axis=1)
    primer_dist_df.columns = ['primer_distance_{}'.format(i) for i in primer_dist_df.columns]
    
    continue_dist_df = pd.DataFrame(continue_distance_counter, index=[sample]).sort_index(axis=1)
    continue_dist_df.columns = ['continue_distance_{}'.format(i) for i in continue_dist_df.columns] 
    
    full_df = full_df.join(primer_dist_df)
    full_df = full_df.join(continue_dist_df)

    full_df.to_csv(report)

    bam_in.close()
    bam_out.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Filter bam file (query name sorted)')
    parser.add_argument('-i', '--input', type=str)
    parser.add_argument('-o','--output', type=str)
    parser.add_argument('-r','--report', type=str)
    parser.add_argument('-s','--sample', type=str)
    
    args = parser.parse_args()

    filter_bam(args.input, args.output, args.report, args.sample)
