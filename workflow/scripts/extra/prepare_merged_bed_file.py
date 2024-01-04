import pandas as pd
import argparse

def read_files(bed_file, insertion_filelist):
    bed_df = pd.read_csv(bed_file, sep='\t', index_col=None, header=None)

    insertion_dfs = [pd.read_csv(ins_file, index_col=0) for ins_file in insertion_filelist]

    return bed_df, insertion_dfs
def merge_matches(name, df):
    ins_name_list = []
    if df.shape[0] == 1:
        series = df.squeeze()
        series[3] = 1
        series.name = name
        ins_name_list.append(series[3])
        return series, ins_name_list
    else:
        res_dict = {}
        res_dict[0] = df[0].unique()[0]
        res_dict[1] = df[1].min()
        res_dict[2] = df[2].max()
        res_dict[3] = len(df[3])
        res_dict[4] = 0
        res_dict[5] = df[5].unique()[0]
        res_dict[6] = df[6]
        series = pd.Series(res_dict)
        series.name = name
        for ins_name in df[3]:
            ins_name_list.append(ins_name)
        return series, ins_name_list

def collapse_bed(bed_df, sample):
    collapsed_bed_df_list = []
    ins_name_to_group_dict_full = {}
    for name, df in bed_df.groupby(6):
        series, ins_name_list = merge_matches(name,df)
        collapsed_bed_df_list.append(series)    
        for ins_name in ins_name_list:
            ins_name_to_group_dict_full[ins_name] = '{}_{}'.format(sample,name)
    collapsed_bed_df = pd.DataFrame(collapsed_bed_df_list)
    return collapsed_bed_df, ins_name_to_group_dict_full

def main():
    parser = argparse.ArgumentParser(description='Analyze merged bed file together with insertion csvs')
    parser.add_argument('--bed', help='Merged bed file')
    parser.add_argument('--insertions', nargs='+', help='List of insertion csvs')
    parser.add_argument('--sample', help='Sample prefix')
    parser.add_argument('--output', help='Output bed file')

    args = parser.parse_args()

    bed_df, insertion_dfs = read_files(args.bed, args.insertions)

    collapsed_bed_df, ins_name_to_group = collapse_bed(bed_df,args.sample)
    print(ins_name_to_group)
    set_of_matches = set()
    for ins_df in insertion_dfs:
        for ins_name in ins_df.index:
            group = ins_name_to_group[ins_name]
            set_of_matches.append([group])
    matches_list = list(set_of_matches)
    collapsed_bed_df_subset = collapsed_bed_df.reindex(matches_list)

    collapsed_bed_df_subset.to_csv(args.output, sep='\t', header=None, index=None)

if __name__ == "__main__":
    main()

