import pysam
import argparse
import matplotlib.pyplot as plt
import mgzip
def quantify_softclip(read):
    cigtup = read.cigartuples
    softclip_len_right = cigtup[-1][1] if cigtup[-1][0] == 4 else 0
    softclip_len_left = cigtup[0][1] if cigtup[0][0] == 4 else 0

    if read.is_reverse:
        return softclip_len_right, softclip_len_left
    else:
        return softclip_len_left, softclip_len_right
def strip_T(seq, qual):
    original_length = len(seq)
    seq = seq.rstrip('T')
    new_length = len(seq)
    R = original_length - new_length
    qual = qual[:new_length]
    seq = seq.lstrip('T')
    final_length = len(seq)
    L = new_length - final_length
    qual = qual[L:]
    return seq, qual, L, R
def format_tuple(tup):
    return b'@'+tup[0].encode()+b'\n'+tup[1].encode()+b'\n+\n'+tup[2].encode()+b'\n'
def extract_transduction_reads(bamfile, file_out_read1):
    f_1 = mgzip.open(file_out_read1, "wb", thread=1,blocksize=2*10**8)
    bam = pysam.AlignmentFile(bamfile,'rb')
    for read in bam.fetch(until_eof=True):
        if read.is_supplementary or read.is_unmapped:
            continue
        softclip_len_start, softclip_len_end = quantify_softclip(read)

        if softclip_len_start >  9:
            forward_sequence = read.get_forward_sequence()
            forward_qual = ''.join([chr(33+i) for i in read.get_forward_qualities()])
            seq, qual, L, R = strip_T(forward_sequence[:softclip_len_start], forward_qual[:softclip_len_start])
            name = '{}-BC:{}-FROM-{}:{}-L:{}-R:{}-Q:{}'.format(read.query_name, read.get_tag('BC'), read.reference_name, read.reference_start, L, R, read.mapping_quality)            
            if len(seq) > 0:
                f_1.write(format_tuple((name, seq, qual)))
        if softclip_len_end >  9:
            forward_sequence = read.get_forward_sequence()
            forward_qual = ''.join([chr(33+i) for i in read.get_forward_qualities()])
            seq, qual, L, R = strip_T(forward_sequence[-softclip_len_end:], forward_qual[-softclip_len_end:])
            name = '{}-BC:{}-TO-{}:{}-L:{}-R:{}-Q:{}'.format(read.query_name, read.get_tag('BC'), read.reference_name, read.reference_start, L, R, read.mapping_quality)            
            if len(seq) > 0:
                f_1.write(format_tuple((name, seq, qual)))
    f_1.close()
    bam.close()
    return None
def main():
    parser = argparse.ArgumentParser(description='Extract polyA-rich reads', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i','--input',metavar='input', type=str, help='Input .bam file')
    parser.add_argument('-o','--output', metavar='output', type=str, help='Output fastq.gz file')
    args = parser.parse_args()

    bamfile = args.input
    file_out_read1 = args.output

    extract_transduction_reads(bamfile, file_out_read1)

if __name__ == '__main__':
    main()
