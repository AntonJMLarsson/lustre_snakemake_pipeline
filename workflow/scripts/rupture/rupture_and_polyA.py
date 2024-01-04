import pickle
import sys
import ruptures as rpt
import numpy as np
import pysam
from collections import Counter

def main():
    loci = sys.argv[1]
    bamfile = sys.argv[2]
    cellfile = sys.argv[3]
    
    cells = set([line.strip() for line in open(cellfile)])
    
    loci_list = loci.split('_')
    if len(loci_list) == 3:
        chrom = loci_list[0]
        start = int(loci_list[1])
        end = int(loci_list[2])
    else:
        chrom = '_'.join(loci_list[:-2])
        start = int(loci_list[-2])
        end = int(loci_list[-1])
    bam = pysam.AlignmentFile(bamfile, 'rb')
    cell_dict = {}
    for read in bam.fetch(chrom,start,end):
        if read.has_tag('LQ'):
            cell = read.get_tag('BC')
            if cell not in cells:
                continue
            a = np.array([ord(c)-33 for c in read.get_tag('LQ')])
            a = 1-10**(-a/10)
            algo = rpt.Binseg(model='l1', jump=1).fit(a)
            result = algo.predict(n_bkps=5)  # Choose the number of change-points
            s = read.get_tag('LS')[8:result[0]]
            if len(s) == 0:
                continue
            if s.count('A')/len(s) < 0.85:
                continue
            if cell not in cell_dict:
                cell_dict[cell] = Counter()
            cell_dict[cell].update({len(s):1})

    with open('polyA_pkls/{}_polyA.pkl'.format(loci), 'wb') as f:
        pickle.dump(cell_dict, f)

if __name__ == "__main__":
    main()
    
