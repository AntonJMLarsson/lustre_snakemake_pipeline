import argparse
import pysam
import numpy as np
ints = {str(i) for i in range(10)}
ops = {'M':0,'I':1, 'D':2, 'N':3, 'S':4, 'H':5,'P':6,'=':7, 'X':8, 'B':9 }
def parse_cigar(cigarstring):
    n = ''
    cigartuple = []
    for c in cigarstring:
        if c in ints:
            n += c
        else:
            cigartuple.append((ops[c], int(n)))
            n = ''
    return cigartuple

def check_read2_softclip(read):
    if read.is_reverse and not read.mate_is_reverse:
        read2_cigar = parse_cigar(read.get_tag('MC'))
        if read2_cigar[0][0] == 4:
            if read2_cigar[0][1] > 23:
                return True
            else:
                return False
        else:
            return False
    elif not read.is_reverse and read.mate_is_reverse:
        read2_cigar = parse_cigar(read.get_tag('MC'))
        if read2_cigar[-1][0] == 4:
            if read2_cigar[-1][1] > 23:
                return True
            else:
                return False
        else:
            return False
    else:
        return False

def check_disconcordant(read, l):
    min_distance = max(l, 15000)
    discond = False
    if read.next_reference_id != read.reference_id:
        discond = True
    elif np.absolute(read.next_reference_start - read.reference_start) > min_distance:
        discond = True
    elif read.is_proper_pair:
        discond = check_read2_softclip(read)
    return discond

def main(in_bamfile, out_bamfile_discon, out_bamfile_concord, l):
    
    bam_in = pysam.AlignmentFile(in_bamfile,'rb')
    
    bam_out_discon = pysam.AlignmentFile(out_bamfile_discon, 'wb', template=bam_in)
    bam_out_concord = pysam.AlignmentFile(out_bamfile_concord, 'wb', template=bam_in)

    for read in bam_in.fetch(until_eof=True):
        if check_disconcordant(read, l):
            bam_out_discon.write(read)
        else:
            bam_out_concord.write(read)

    return True

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Remove concordant reads.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i','--input',metavar='input', type=str, help='Input .bam file')
    parser.add_argument('-o1','--output_discon', metavar='output', type=str, help='Output .bam file')
    parser.add_argument('-o2','--output_concord', metavar='output', type=str, help='Output .bam file')
    parser.add_argument('-l','--length', metavar='length', type=int, default=15000, help='Distance on the same chromosome to be considered disconcordant')

    args = parser.parse_args()
    
    main(args.input, args.output_discon, args.output_concord, args.length)

