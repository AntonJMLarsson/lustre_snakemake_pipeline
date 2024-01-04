import sys
import pysam
def main():
    bam_in = pysam.AlignmentFile(sys.argv[1], 'rb')
    bam_out = pysam.AlignmentFile(sys.argv[2], 'wb', template=bam_in)

    for read in bam_in.fetch(until_eof=True):
        readidlist = read.query_name.split("-")
        read.query_name = readidlist[0]
        read.set_tag(tag = "BC", value_type = "Z", value = readidlist[1].replace("BC:",""))
        bam_out.write(read)
    bam_in.close()
    bam_out.close()

if __name__ == "__main__":
    main()
