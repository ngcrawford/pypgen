#!/usr/bin/env python
# encoding: utf-8

import os
import re
import sys
import datetime
import argparse
from VCF import *
from helpers import *
import multiprocessing
from collections import defaultdict

"""
python vcf2sweepfinder.py \
-i test_data/butterfly.vcf.gz \
-o bayesscan.snps \
-p cydno:c511,c512,c513,c514,c515,c563,c614,c630,c639,c640 \
outgroup:h665,i02-210 \
melpo:m523,m524,m525,m589,m675,m676,m682,m683,m687,m689 \
pachi:p516,p517,p518,p519,p520,p591,p596,p690,p694,p696 \
-c 2 \
-r Chr01:1-10001
"""

def vcf_line_2_sweepfinder_line(vcf_line, populations, outgroup_call):

    """
    The snp file should be a tab-delimited file with column headers, and
one row per SNP.  One column header should be "x" (the frequency of the
SNP), another should be "n" (the sample size, must be greater than x), and
another should be "position" (the chromosomal location of the SNP).  
Optionally, a fourth column named "folded" can be added.  If it is present,
than a value of one indicates that the SNP is folded (there is no distinction
between ancestral/derived states), and 0 means unfolded.  If the folded
column is not present, all SNPs are assumed to be unfolded.
Column names do not actually contain quotes.  A sample input file might look
something like:

position        x       n   folded
37.000000       10      46  0
145.000000      3       47  0
277.000000      1       47  1
385.000000      37      43  1
469.000000      2       45  0
585.000000      1       44  0
733.000000      10      45  0

    """
    
    # Process Folding
    if outgroup_call == None: 
        folded = 0
    else: 
        folded = 1


    print calc_allele_counts(populations, vcf_line)



def main():
    # get args. 
    args = default_args()
    args = args.parse_args()
    
    # TODO: 
    # test that pysam is installed.
    # bgzip check. MDSum?

    # Get slice indicies
    # 1. read file and get chrm sizes
    # 2. process chrm sizes and return as
    #    slices and a zipped list (chrm, (start, stop))
    slice_indicies = get_slice_indicies(args.input, 
        args.regions, args.window_size, args.regions_to_skip)
    

    starting_time = datetime.datetime.now()
    previous_update_time = datetime.datetime.now()
    
    # Calculate the total size of the dataset
    chrm_lengths = process_header(args.input)
    total_bp_in_dataset = sum(chrm_lengths.values())

    # Get information about samples from the header.
    # this becomes the precursor to the VCF row
    empty_vcf_line = make_empty_vcf_ordered_dict(args.input)


    # Convert populations input into a dict of pops where
    # values are lists of samples
    populations = parse_populations_list(args.populations)

    order = [] # store order of samples.
    bp_processed = 0

    p = multiprocessing.Pool(processes=int(args.cores),maxtasksperchild=1000)

    fstat_input_iterator = generate_fstats_from_vcf_slices(slice_indicies, populations, empty_vcf_line, args)
    
    for fs in fstat_input_iterator:
        
        tabix_slice, chrm, start, stop, populations, header, min_samples = fs
        for count, line in enumerate(tabix_slice):
            vcf_line = parse_vcf_line(line, header)

            if vcf_line["FILTER"] != "PASS": continue

            outgroup_call = process_outgroup(vcf_line, populations)
            print vcf_line_2_sweepfinder_line(vcf_line, populations, outgroup_call)

        break

if __name__ == '__main__':
    main()




