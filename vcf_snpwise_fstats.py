#!/usr/bin/env python
# encoding: utf-8

"""
python sliding_window.py \
-i test_data/butterfly.vcf.gz \
-o bayesscan.snps \
-p cydno:c511,c512,c513,c514,c515,c563,c614,c630,c639,c640 \
outgroups:h665,i02-210 \
melpo:m523,m524,m525,m589,m675,m676,m682,m683,m687,m689 \
pachi:p516,p517,p518,p519,p520,p591,p596,p690,p694,p696 \
-c 2 \
-r Chr01:1-10001
"""


import os
import re
import sys
import gzip
import datetime
import argparse
from VCF import *
from helpers import *
import multiprocessing


def process_header(tabix_file):

    chrm_lenghts_dict = {}
    tabix_file = pysam.Tabixfile(tabix_file)
    
    for line in tabix_file.header:
        
        if line.startswith("##contig") == True: 
            chrm, length = re.split(r"<ID=|,length=", line)[1:]
            length =  int(length.strip(">"))
            chrm_lenghts_dict[chrm] = length

    return chrm_lenghts_dict


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
    slice_indicies = get_slice_indicies(args.input, args.regions, args.window_size, args.regions_to_skip)
    starting_time = datetime.datetime.now()
    previous_update_time = datetime.datetime.now()
    bp_processed = 0
    
    # Calculate the total size of the dataset
    chrm_lengths = process_header(args.input)
    total_bp_in_dataset = sum(chrm_lengths.values())

    # Get information about samples from the header.
    # this becomes the precursor to the VCF row
    empty_vcf_line = make_empty_vcf_ordered_dict(args.input)

    # Convert populations input into a dict of pops where
    # values are lists of samples
    populations = parse_populations_list(args.populations)

    fstat_order = [] # store order of paired samples.
    pop_size_order = []
    vcf_count = 0


    for count, line in enumerate(open_vcf(args)):
        
        if line.startswith('#') == True: continue
        vcf_count += 1

        vcf_line = parse_vcf_line(line, empty_vcf_line)
        allele_counts = calc_allele_counts(populations, vcf_line)
        fstats =  calc_fstats(allele_counts)
        pop_size_stats = get_population_sizes(vcf_line, populations)

        f_stats, fstat_order = f_statistics_2_sorted_list(fstats, order=fstat_order)
        pop_size_stats, pop_size_order = pop_size_statistics_2_sorted_list(pop_size_stats, order=pop_size_order)

        chrm = vcf_line['CHROM']
        pos = vcf_line['POS']

        if vcf_count == 1:
            args.output.write(','.join(['chrm','pos','snp_count', 'total_depth_mean', 'total_depth_stdev'] \
                       + map(str, pop_size_order + fstat_order) + "\n")
            args.output.write(','.join(map(str, [chrm, pos] + pop_size_stats + f_stats)) + "\n")

        else:
            args.output.write(','.join(map(str, [chrm, pos] + pop_size_stats + f_stats)) + "\n")


        #previous_update_time = progress_meter(previous_update_time, chrm, stop, bp_processed, total_bp_in_dataset)


if __name__ == '__main__':
    main()
