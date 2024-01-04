import pysam
import networkx as nx
import pandas as pd
from collections import Counter
import argparse

def find_overlap(val, intervals):
    for idx, interval in intervals.items():
        if val >= interval[0] and val <= interval[1]:
            return idx
    return None

def have_bidirectional_relationship(G, node1, node2):
    return G.has_edge(node1, node2) and G.has_edge(node2, node1)


def main():
    parser = argparse.ArgumentParser(description='Count KR cell detection', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input',metavar='BAM', type=str, help='Input .bam file')
    parser.add_argument('-o', '--output',metavar='CSV', type=str, help='Output csv file')
    parser.add_argument('--KR', metavar='KR', type=str, help='Input KR insertion bed file')
    parser.add_argument('--KNR', metavar='KNR', type=str, help='Input KNR insertion bed file')
    parser.add_argument('--UNK', metavar='UNK', type=str, help='Input UNK insertion bed file')
    args = parser.parse_args()

    bamfile = args.input
    outfile = args.output
    KNR_file = args.KNR
    KR_file = args.KR
    UNK_file = args.UNK

    df_KNR = pd.read_csv(KNR_file, sep='\t', header=None)
    df_KNR[1] = df_KNR[1]-1000
    df_KNR[2] = df_KNR[2]+1000
    df_KNR['type'] = 'KNR'
    df_KR = pd.read_csv(KR_file, sep='\t', header=None)
    df_KR[1] = df_KR[1]-1000
    df_KR[2] = df_KR[2]+1000
    df_KR['type'] = 'KR'

    bam = pysam.AlignmentFile(bamfile, 'rb')
    contigs = set([d['SN'] for d in bam.header['SQ']])

    df_UNK = pd.read_csv(UNK_file, index_col=0)
    df_UNK.columns = [0,1,2]
    df_UNK['type'] = 'UNK'
    df_all = df_KNR.append(df_KR).append(df_UNK).reset_index()

    L1_intervals = {}
    for idx, row in df_all.iterrows():
        if row[0] not in L1_intervals:
            L1_intervals[row[0]] = {}
        L1_intervals[row[0]][idx] = [row[1], row[2]]
    
    direction_dict = {}
    G = nx.DiGraph()
    for idx, row in df_all.iterrows():
        if row[0] not in contigs:
            continue
        idx_counter = Counter()
        direction_dict_bla = Counter()
        for read in bam.fetch(row[0], max(0,row[1]), row[2]):
            if read.mapping_quality < 25:
                continue
            query_list = read.query_name.split('-')
            coords_list = query_list[3].split(':')
            chrom_other = coords_list[0]
            start_other = int(coords_list[1])
            if int(query_list[6].split(':')[1]) < 25:
                continue
            if chrom_other not in L1_intervals:
                continue
            source = find_overlap(start_other, L1_intervals[chrom_other])
            if source is None:
                continue
            idx_counter.update({source:1})
            if source not in direction_dict_bla:
                direction_dict_bla[source] = Counter()
            direction_dict_bla[source].update({query_list[2]:1})
        sources = {source for source, n in idx_counter.items() if n > 10}
        if len(sources) > 0:
            direction_dict[idx] = {}
            for source in sources:
                direction_dict[idx][source] = direction_dict_bla[source].most_common()[0][0]
                G.add_edge(idx, source)
    biconnections = set()
    for u, v in G.edges():
        if u > v:  # Avoid duplicates, such as (1, 2) and (2, 1)
            v, u = u, v
        if have_bidirectional_relationship(G, u, v):
            biconnections.add((u, v))
    

    result = []
    for t in biconnections:
        d_1 = direction_dict[t[0]][t[1]]
        d_2 = direction_dict[t[1]][t[0]]
        if d_1 == d_2:
            continue
        ins_1 = df_all.loc[t[0]]
        ins_2 = df_all.loc[t[1]]
        if t[0] > 16855:
            print(d_1,t)
        if t[1] > 16855:
            print(d_1,t)
        if d_1 == 'TO':
            result.append([ins_1[0], ins_1[1], ins_1[2],ins_1['type'],ins_2[0], ins_2[1], ins_2[2], ins_2['type']])
        else:
            result.append([ins_2[0], ins_2[1], ins_2[2],ins_2['type'], ins_1[0], ins_1[1], ins_1[2], ins_1['type']])
    
    df_result = pd.DataFrame(result)
    df_result.columns = ['chr1','start1','end1','source_type','chr2','start2','end2','target_type']
    df_result.to_csv(outfile, index=None)

if __name__ == "__main__":
    main()