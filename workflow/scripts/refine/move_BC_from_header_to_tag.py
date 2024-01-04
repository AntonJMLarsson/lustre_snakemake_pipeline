import pysam
import argparse
def move_from_header_to_tag(bamfile_in, bamfile_out, tag):
    bam_in = pysam.AlignmentFile(bamfile_in, 'rb')
    bam_out = pysam.AlignmentFile(bamfile_out, 'wb', template=bam_in)

    BC = bam_in.header['RG'][0][tag]
    for read in bam_in.fetch(until_eof=True):
        read.set_tag('BC', BC)
        bam_out.write(read)

    bam_in.close()
    bam_out.close()

def main():
    parser = argparse.ArgumentParser(description='Move cell barcode from header into tag', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i','--input',metavar='input', type=str, help='Input .bam file')
    parser.add_argument('-o','--output', metavar='output', type=str, help='Output .bam file')
    parser.add_argument('--tag', type=str, default='BC', help='Name of tag in header')

    args = parser.parse_args()

    bamfile_in = args.input
    bamfile_out = args.output
    tag = args.tag

    move_from_header_to_tag(bamfile_in, bamfile_out, tag)

if __name__ == '__main__':
    main()
