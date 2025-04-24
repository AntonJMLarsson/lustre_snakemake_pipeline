import pysam
import argparse
def convert_custom_coordinates(read):
    ref_name = read.reference_name
    if ref_name is None:
        return (None, -2, 0)
    ref_name_list = ref_name.split('_')

    if len(ref_name_list) == 3:
        chrom = ref_name_list[0]
    else:
        chrom = '_'.join(ref_name_list[:-2])
    start = int(ref_name_list[-2])+read.reference_start-1000
    return (chrom, start, read.mapping_quality)
def transform_reads(bamfile_in, header_bam, bamfile_out):
    bam_in = pysam.AlignmentFile(bamfile_in,'rb')
    bam_h = pysam.AlignmentFile(header_bam,'rb')
    bam_out = pysam.AlignmentFile(bamfile_out, 'wb',template=bam_h)
    cell = bamfile_in.split('/')[-1].split('_')[0]

    for read in bam_in.fetch(until_eof=True):
        read_dict = read.to_dict()
        tup = convert_custom_coordinates(read)
        read_dict['ref_name'] = str(tup[0])
        read_dict['ref_pos'] = str(tup[1]+1)
        if tup[0] is None:
            print('Unmapped')
        read = pysam.AlignedSegment.from_dict(read_dict, bam_out.header)
        read.set_tag('BC', cell)
        bam_out.write(read)

    bam_in.close()
    bam_out.close()

def main():
    parser = argparse.ArgumentParser(description='Move cell barcode from header into tag', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i','--input',metavar='input', type=str, help='Input .bam file')
    parser.add_argument('-head','--header',metavar='input', type=str, help='Header .bam file')
    parser.add_argument('-o','--output', metavar='output', type=str, help='Output .bam file')

    args = parser.parse_args()

    bamfile_in = args.input
    header_bam = args.header
    bamfile_out = args.output

    transform_reads(bamfile_in, header_bam, bamfile_out)

if __name__ == '__main__':
    main()

