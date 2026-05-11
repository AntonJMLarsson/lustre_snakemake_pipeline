import pysam
import argparse
import matplotlib.pyplot as plt
import mgzip
def quantify_A_content_and_softclip(read):
    aligned_sequence = read.query_alignment_sequence
    cigtup = read.cigartuples
    if read.is_reverse:
        fraction_As = aligned_sequence.count('A')/len(aligned_sequence)
        softclip_len = cigtup[-1][1] if cigtup[-1][0] == 4 else 0
    else:
        fraction_As = aligned_sequence.count('T')/len(aligned_sequence)
        softclip_len = cigtup[0][1] if cigtup[0][0] == 4 else 0
    return fraction_As, softclip_len
def format_tuple(tup):
    return b'@'+tup[0].encode()+b'\n'+tup[1].encode()+b'\n+\n'+tup[2].encode()+b'\n'
def extract_polyA_rich_reads(bamfile, file_out_read1):
    f_1 = mgzip.open(file_out_read1, "wb", thread=1,blocksize=2*10**8)
    bam = pysam.AlignmentFile(bamfile,'rb')
    fraction_As_list = []
    softclip_len_list = []
    f_1.write(format_tuple(('dummy', 'AAA', 'AAA')))
    for read in bam.fetch(until_eof=True):
        if read.is_supplementary or read.is_unmapped:
            continue
        fraction_As, softclip_len = quantify_A_content_and_softclip(read)

        if fraction_As > 0.9 and softclip_len > 0:
            forward_sequence = read.get_forward_sequence()
            forward_qual = ''.join([chr(33+i) for i in read.get_forward_qualities()])
            padding = min(10, len(forward_sequence)-softclip_len)
            f_1.write(format_tuple((read.query_name, forward_sequence[:softclip_len+padding], forward_qual[:softclip_len+padding])))
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

    extract_polyA_rich_reads(bamfile, file_out_read1)

if __name__ == '__main__':
    main()
