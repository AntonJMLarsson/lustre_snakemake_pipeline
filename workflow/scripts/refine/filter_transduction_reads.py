import pysam
import argparse
def filter_transduction_reads(bamfile_in, bamfile_out):
    bam_in = pysam.AlignmentFile(bamfile_in, 'rb')
    bam_out = pysam.AlignmentFile(bamfile_out, 'wb', template=bam_in)

    for read in bam_in.fetch(until_eof=True):
        if read.mapping_quality == 0:
            continue
        query_list = read.query_name.split('-')
        if query_list[6] == 'Q:0':
            continue
        coords_list = query_list[3].split(':')
        chrom_other = coords_list[0]
        start_other = int(coords_list[1])
        if chrom_other == read.reference_name and abs(start_other - read.reference_start) < 10000:
            continue
        bam_out.write(read)
    bam_in.close()
    bam_out.close()

def main():
    parser = argparse.ArgumentParser(description='Filter transduction reads', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i','--input',metavar='input', type=str, help='Input .bam file')
    parser.add_argument('-o','--output', metavar='output', type=str, help='Output .bam file')
    args = parser.parse_args()

    bamfile_in = args.input
    bamfile_out = args.output

    filter_transduction_reads(bamfile_in, bamfile_out)

if __name__ == '__main__':
    main()
