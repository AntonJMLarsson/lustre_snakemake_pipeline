import pysam
import argparse
def modify_read(read1, seq, qual, cs, flag):
    read1.set_tag('LS',seq[23:])
    qual = ''.join([chr(33+i) for i in qual])
    read1.set_tag('LQ',qual[23:])
    read1.set_tag('CS',cs)
    read1.set_tag('FL', flag)
    return read1

def split_file(bamfile, outfile_r1, outfile_r2):
    bam = pysam.AlignmentFile(bamfile, 'rb')
    bam_r1 = pysam.AlignmentFile(outfile_r1, 'wb', template=bam)
    bam_r2 = pysam.AlignmentFile(outfile_r2, 'wb', template=bam)

    for read in bam:
        if read.is_read1 and read.is_paired and not read.is_supplementary:
            read2 = next(bam)
            if read2.is_supplementary:
                while read2.is_supplementary:
                    read2 = next(bam)
            assert read.query_name == read2.query_name
            seq = read2.get_forward_sequence()
            qual = read2.get_forward_qualities()
            cs = read2.cigarstring
            flag = read2.flag
            read = modify_read(read, seq, qual, cs, flag)
            bam_r1.write(read)
    bam_r1.close()
    bam_r2.close()
    return None

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Split bam file into read1 and read2')
    parser.add_argument('-i', '--input', type=str)
    parser.add_argument('-r1','--out-r1', type=str)
    parser.add_argument('-r2','--out-r2', type=str)

    args = parser.parse_args()

    split_file(args.input, args.out_r1, args.out_r2)


    
