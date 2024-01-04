import pandas as pd
import argparse

def write_bed_file(file_list, bed_file):
    df_list = []
    for filename in file_list:
        df = pd.read_csv(filename, index_col=0)
        df_list.append(df)
    full_df = pd.concat(df_list)
    full_df[['chrom', 'start','end']].to_csv(bed_file, sep='\t', header=None, index=False)

def main():
    parser = argparse.ArgumentParser(description='Make bed file', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i','--input', nargs='+', metavar='list', required=True, help='List of input files')
    parser.add_argument('-o','--output', metavar='bed', type=str, required=True, help='Output bed file')
    args = parser.parse_args()
    
    file_list = args.input
    bed_file = args.output

    print("Writing genomic coordinates from {} to {}".format(file_list, bed_file))

    write_bed_file(file_list, bed_file)

if __name__ == "__main__":
    main()
