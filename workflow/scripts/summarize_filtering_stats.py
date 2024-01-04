import glob
import pandas as pd
import numpy as np
import argparse

def summarize_stats(folder, samplesheet_file, outfile):
    files = glob.glob('{}/*.csv'.format(folder))
    dfs = [pd.read_csv(file, index_col=0) for file in files]
    full_df = pd.concat(dfs)

    if samplesheet_file != 'None':
        samplesheet = pd.read_csv(samplesheet_file, index_col=0)
        full_df = full_df.join(samplesheet)

    full_df['total_n'] = full_df['n_pass']+full_df['n_fail']
    full_df['fraction_pass'] = full_df['n_pass']/full_df['total_n']
    
    full_df.to_csv(outfile)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Concatenate stats output from many individual files', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--folder', type=str)
    parser.add_argument('-s','--samplesheet', type=str)
    parser.add_argument('-o','--outfile',type=str)

    args = parser.parse_args()
    
    summarize_stats(args.folder, args.samplesheet, args.outfile)
