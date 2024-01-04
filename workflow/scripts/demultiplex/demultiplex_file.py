import mgzip
from fastqandfurious import fastqandfurious
import argparse
import pysam
import pandas as pd
import Levenshtein
from contextlib import contextmanager
def modify_read(read1, read2, bc):
    seq, qual, flag = read2.get_forward_sequence(), read2.get_forward_qualities(), read2.flag
    read1.set_tag('BC',bc)
    read1.set_tag('LS',seq)
    qual = ''.join([chr(33+i) for i in qual])
    read1.set_tag('LQ',qual)
    read1.set_tag('FL', flag)
    return read1
def get_read_name(r2, illumina):
    if illumina:
        name = r2[0].decode('utf-8').split(' ')[0]
    else:
        name = r2[0][:-2].decode('utf-8')
    return name

import fastqandfurious._fastqandfurious as _fqf
from fastqandfurious.fastqandfurious import entryfunc

@contextmanager
def multi_file_manager(files, bamfile, mode='rt'):
	"""
	Open multiple files and make sure they all get closed.
	"""
	temp = pysam.AlignmentFile(bamfile, "rb")
	files = [pysam.AlignmentFile(f, "wb", template = temp) for f in files]
	temp.close()
	yield files
	for f in files:
		f.close()

def add_BC(r2fq_file, metadata_file, corrected_bc_file, bam_infile, outfolder, illumina, number):
    def custom_entryfunc(buf, pos, globaloffset):
        header = buf[(pos[0]+1):pos[1]]
        sequence = buf[pos[2]:pos[3]]
        return (header, sequence)

    pin = 10
    bufsize = 20000
    r2fq = r2fq_file
    r2_infile = mgzip.open(r2fq, 'rb', thread = pin)
    r2_it = fastqandfurious.readfastq_iter(r2_infile, bufsize, entryfunc=custom_entryfunc)
    r2 = next(r2_it)

    rname = get_read_name(r2, illumina)
    bc = r2[1][-20:].decode('utf-8')


    metadata_df = pd.read_csv(metadata_file, index_col = 0)

    corrected_bcs_df = pd.read_csv(corrected_bc_file, index_col = 0)

    cell_bcs = set(metadata_df[metadata_df['Modality'] == 'DNA'].index)

    corrected_bcs_df = corrected_bcs_df.reindex(corrected_bcs_df.index[corrected_bcs_df.apply(lambda row: row['trueBC'] in cell_bcs, axis=1)])

    falseBC_to_trueBC = {}
    for falseBC, row in corrected_bcs_df.iterrows():
        falseBC_to_trueBC[falseBC] = row['trueBC']

    bam_qname_sorted = pysam.AlignmentFile(bam_infile, 'r', ignore_truncation=True)


    aligned_r1 = {}
    aligned_r2 = {}
    read_bcs = {}
    aligned_supplementary_list = {}
    read_bcs[rname] = bc
    print(rname, bc)
    print(bc in metadata_df.index)
    t = True
    i = 0

    files = [outfolder + "/" + bc1 + '.' + str(number) + ".bam" for bc1 in cell_bcs]
    bc_dict = {bc1: i  for i,bc1 in enumerate(cell_bcs)}

    def yield_reads(bam_qname_sorted):
        try:
            for aligned_read in bam_qname_sorted.fetch(until_eof=True):
                yield aligned_read
        except OSError:
            return None
        return None
    with multi_file_manager(files, bam_infile, mode='rt') as fopen:
        for aligned_read in yield_reads(bam_qname_sorted):
            if aligned_read is None:
                break
            if i == 0:
                prev_qname = aligned_read.query_name
                while aligned_read.query_name != rname:
                    r2 = next(r2_it, None)
                    rname = get_read_name(r2, illumina)
                bc = r2[1][-20:].decode('utf-8')
                read_bcs[rname] = bc
                i = 1
            if t:
                print(aligned_read.query_name)
                t = False
    
            if aligned_read.is_read1 and not aligned_read.is_supplementary and aligned_read.is_paired:
                aligned_r1[aligned_read.query_name] = aligned_read
            elif aligned_read.is_read2 and not aligned_read.is_supplementary and aligned_read.is_paired:
                aligned_r2[aligned_read.query_name] = aligned_read
            elif aligned_read.is_supplementary:
                if aligned_read.query_name in aligned_supplementary_list:
                    aligned_supplementary_list[aligned_read.query_name].append(aligned_read)
                else:
                    aligned_supplementary_list[aligned_read.query_name] = [aligned_read]
    
            while aligned_read.query_name not in read_bcs:
                r2 = next(r2_it, None)
                read_bcs[get_read_name(r2, illumina)] = r2[1][-20:].decode('utf-8')
            if prev_qname != aligned_read.query_name:
                BC_tag = 'bla'
                if read_bcs[prev_qname] in falseBC_to_trueBC:
                    BC_tag = falseBC_to_trueBC[read_bcs[prev_qname]]
                else:
                    BC_tag = read_bcs[prev_qname]
                
                
                
                if BC_tag in cell_bcs:
                    f_BC = fopen[bc_dict[BC_tag]]
                    f_BC.write(modify_read(aligned_r1[prev_qname], aligned_r2[prev_qname], BC_tag))
                    f_BC.write(modify_read(aligned_r2[prev_qname], aligned_r1[prev_qname], BC_tag))
                    if prev_qname in aligned_supplementary_list:
                        for supp_read in aligned_supplementary_list[prev_qname]:
                            supp_read.set_tag('BC', BC_tag)
                            f_BC.write(supp_read)
                del aligned_r1[prev_qname]
                del aligned_r2[prev_qname]
                if prev_qname in aligned_supplementary_list:
                    del aligned_supplementary_list[prev_qname]
                del read_bcs[prev_qname]
            prev_qname = aligned_read.query_name
    return None

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Merge mapped BWA file with an unmapped demultiplexed bam file')
    parser.add_argument('-m','--mapped', metavar='mapped',  type=str, help='Mapped BWA file, query name sorted')
    parser.add_argument('-f','--fastq', metavar='fastq', type=str, help='fastq file with cell barcodes')
    parser.add_argument('-o','--outfolder', metavar='outfolder', type=str, help='Output directory')
    parser.add_argument('-b', '--barcodes', metavar='barcodes', type=str, help='Barcode file from zUMIs with correction')
    parser.add_argument('-n', '--number', metavar='number', type=str, help='Number to append')
    parser.add_argument('-meta','--metadata', metavar='metadata', type=str, help='Metadata of cell barcodes')
    parser.add_argument('--illumina', action='store_true', help='Indicate if experiment is sequenced with Illumina')
    args = parser.parse_args()

    fastq_file = args.fastq
    mapped_bam_file = args.mapped
    outfolder = args.outfolder
    corrected_bc_file = args.barcodes
    number = args.number
    metadata_file = args.metadata
    illumina = args.illumina

    print('Adding barcodes from {} with correction from {} to bam file {} for cells present in {}. Output written to: {}'.format(fastq_file, corrected_bc_file, mapped_bam_file, metadata_file, outfolder))

    add_BC(fastq_file, metadata_file, corrected_bc_file, mapped_bam_file, outfolder, illumina, number)




