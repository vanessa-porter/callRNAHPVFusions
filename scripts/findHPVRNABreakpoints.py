#!/usr/bin/env python
"""detectHpvRnaIntegration.py: Detects """

import yaml
import argparse
import os
import sys
import pandas as pd
from collections import Counter

#file = "output/TCGA-C5-A8ZZ/star/hpvRNAChimeric.out.filtered.junction"
#out = "output/HTMCP-03-06-02108/hpvRNAbreakpoints.txt"

# Check the number of command line arguments
if not len(sys.argv)==3:
    print("\nError:\tincorrect number of command-line arguments")
    print("Syntax:\tlookForSVs.py [Input Junction] [OutputDir]\n")
    sys.exit()

# File input
#with open(file) as f:
#    f = [l.strip().split("\t") for l in f.readlines()]

with open(sys.argv[1]) as f:
    f = [l.strip().split("\t") for l in f.readlines()]

### MAKE A VCF FILE CONTAINING ONLY THE ASSOCIATED SVs
chromosomes = ["chr" + str(i) if i != 'X' else "chrX" for i in range(1, 22)]
bp_list = []
for strLine in f:
    chr1 = strLine[0]
    chr2 = strLine[3]
    if ((chr1.__contains__('HPV') and chr2 in chromosomes) or (chr2.__contains__('HPV') and chr1 in chromosomes)):
        bp1 = strLine[0] + ":" + strLine[1] + ":" + strLine[2]
        bp2 = strLine[3] + ":" + strLine[4] + ":" + strLine[5]
        if (bp1.__contains__('HPV')):
            hpvbp = bp1
            humanbp = bp2
        else:
            hpvbp = bp2
            humanbp = bp1
        bp = humanbp + "^" + hpvbp     
        bp_list.append(bp)

# count the replicate reads with the same breakpoint
counts = Counter(bp_list)

# filter for breakpoints found in > 10 reads
expr_bp = {k: v for k, v in counts.items() if v > 10}

# create data frames for text file and for bed file
bp_df = pd.DataFrame(expr_bp.items(), columns=['breakpair', 'num.of.reads'])
bp_df[['human.bp', 'hpv.bp']] = bp_df['breakpair'].str.split('^', n=1, expand=True)
bp_df[['chr', 'pos', 'strand']] = bp_df['human.bp'].str.split(':', expand=True)
bp_df[['HPVchr', 'HPVpos', 'HPVstrand']] = bp_df['hpv.bp'].str.split(':', expand=True)
bpdf = bp_df[['chr', 'pos', 'strand', 'HPVchr', 'HPVpos', 'HPVstrand', 'num.of.reads']]

# make a bed score using the num of reads (% of fusion reads)
bp_df['num.of.reads'] = bp_df['num.of.reads'].astype(int)
total_reads = bp_df['num.of.reads'].sum()
bp_df['score'] = (bp_df['num.of.reads'] / total_reads) * 1000

#make the bp_bed dataframe
bp_bed = bp_df[['chr', 'pos','strand', 'breakpair', 'score','num.of.reads']]
bp_bed['pos'] = bp_bed['pos'].astype(int)
bp_bed['pos2'] = bp_bed['pos'] + 1
bp_bed = bp_bed[['chr', 'pos', 'pos2', 'breakpair', 'score','strand','num.of.reads']]

# save the dataframes
#out1 = sys.argv[2] + '/hpvRNAbreakpoints.txt'
out2 = sys.argv[2] + '/hpvRNAbreakpoints.bed'
#bpdf.to_csv(path_or_buf = out1, sep = "\t", header = True, index = False)
bp_bed.to_csv(path_or_buf = out2, sep = "\t", header = False, index = False)