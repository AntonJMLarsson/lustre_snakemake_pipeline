from pyfaidx import Fasta
import pandas as pd
import argparse
def get_seq(contig, start, end, fasta_ref):
    string = fasta_ref[contig][start-1000:end+1000].seq.upper()
    return string

def make_fasta(bed_file, fasta_file, out_file):

    fasta_ref = Fasta(fasta_file)
    df = pd.read_csv(bed_file, sep='\t', comment='t', header=None, index_col=None).dropna()
    with open(out_file, 'w') as f_ref:
        for i, row in df.iterrows():
            if row[0] not in fasta_ref:
                continue
            ref_seq = get_seq(row[0], int(row[1]), int(row[2]), fasta_ref)
            f_ref.write('>{}_{}_{}\n'.format(row[0], int(row[1]), int(row[2])))
            f_ref.write('{}\n'.format(ref_seq))
    return True

def main():
    parser = argparse.ArgumentParser(description='Write fasta from bed file', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i','--input',metavar='input', type=str, help='Input .bed file')
    parser.add_argument('-f','--fasta', metavar='fasta', type=str, help='Reference fasta file')
    parser.add_argument('-o','--output', metavar='output', type=str, help='Output .fa file')

    args = parser.parse_args()

    bed_file = args.input
    fasta_file = args.fasta
    out_file = args.output

    make_fasta(bed_file, fasta_file, out_file)

if __name__ == "__main__":
    main()
