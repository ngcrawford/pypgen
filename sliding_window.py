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
import datetime
import argparse
from VCF import *
from helpers import *
import multiprocessing

def generate_fstats_from_vcf_slices(slice_indicies, populations, header, args):

    for count, si in enumerate(slice_indicies):
        chrm, start, stop = si

        yield [slice_vcf(args.input, chrm, start, stop), 
               chrm, start, stop, populations, header, 
               args.min_samples]

        # if count > 1: break

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
    
    p = multiprocessing.Pool(processes=int(args.cores), maxtasksperchild=1000)

    fstat_input_iterator = generate_fstats_from_vcf_slices(slice_indicies, populations, empty_vcf_line, args)
    for count, result in enumerate(p.imap(calc_slice_stats, fstat_input_iterator)):
        
        if result == None: continue
        
        chrm_start_stop, pop_size_statistics, fstats = result
        
        # TO DO: Figure out why some samples have no data (BUG?!)
        if not pop_size_statistics: continue
        if not fstats: continue

        f_stats, fstat_order = multilocus_f_statistics_2_sorted_list(fstats, order=fstat_order)
        pop_size_stats, pop_size_order = pop_size_statistics_2_sorted_list(pop_size_statistics, order=pop_size_order)

        if count == 0:
            args.output.write(','.join(['chrm','start','stop','snp_count', 'total_depth_mean', 'total_depth_stdev'] \
                       + map(str, pop_size_order) + map(str, fstat_order)) + "\n")

            args.output.write(','.join(map(str, chrm_start_stop) + map(str, pop_size_stats) + map(str, f_stats)) + "\n")
        else:
            args.output.write(','.join(map(str, chrm_start_stop) + map(str, pop_size_stats) + map(str, f_stats)) + "\n")

        chrm, start, stop  = chrm_start_stop[0:3]
        bp_processed += args.window_size

        #previous_update_time = progress_meter(previous_update_time, chrm, stop, bp_processed, total_bp_in_dataset)


if __name__ == '__main__':
    main()