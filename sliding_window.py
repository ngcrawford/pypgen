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

def get_args():
    """Parse sys.argv"""
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input', 
                        required=True, 
                        type=str,
                        help='Path to VCF file.')
    
    parser.add_argument('-o','--output',
                        help='Path to output csv file. \
                              If path is not set defaults to STDOUT.')
    
    parser.add_argument('-c','--cores', 
                        required=True, 
                        type=int,
                        help='Number of cores to use.')

    parser.add_argument('-r','--regions', 
                        required=False,
                        # action=FooAction, # TODO: fix this!
                        nargs='+',
                        help="Define a chromosomal region. \
                              A region can be presented, for example, in the following \
                              format: ‘chr2’ (the whole chr2), ‘chr2:1000000’ (region \
                              starting from 1,000,000bp) or ‘chr2:1,000,000-2,000,000’ \
                              (region between 1,000,000 and 2,000,000bp including the end \
                              points). The coordinate is 1-based.' [Same format as \
                              SAMTOOLs/GATK, example text cribbed from SAMTOOLs]")

    parser.add_argument('-f', '--filter-string',
                        require=False,
                        default="'FILTER' == 'PASS'")

    parser.add_argument('--regions-to-skip',
                        default=[],
                        required=False,
                        nargs='+',
                        help='Define a chromosomal region(s) to skip.')
    
    parser.add_argument('-p','--populations', 
                        nargs='+',
                        help='Names of populations and samples. \
                              The format is: "PopName:sample1,sample2,sample3,etc..."')

    parser.add_argument('-w','--window-size',
                        default=5000, 
                        type=int,
                        help='Size of the window in which to \
                              calculate pairwise F-staticstics')

    parser.add_argument("-m",'--min-samples',
                        type=int,
                        default=5,
                        help="Minimum number of samples per population.")

    return parser.parse_args()

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
    args = get_args()
    
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

    fout = open(args.output,'w')    

    p = multiprocessing.Pool(processes=int(args.cores),maxtasksperchild=1000)

    fstat_input_iterator = generate_fstats_from_vcf_slices(slice_indicies, populations, empty_vcf_line, args)
    for count, result in enumerate(p.imap(calc_slice_stats, fstat_input_iterator)):
        
        if result == None: continue
        
        chrm_start_stop, result = result
        stats, order = multilocus_f_statistics_2_sorted_list(result, order=order)

        if count == 0:
            fout.write(','.join(['chrm','start','stop','snp_count', 'total_depth_mean', 'total_depth_stdev'] +  map(str, order)) + "\n")
            fout.write(','.join(map(str, chrm_start_stop) + map(str, stats)) + "\n")
        else:
            fout.write(','.join(map(str, chrm_start_stop) + map(str, stats)) + "\n")

        chrm, start, stop  = chrm_start_stop[0:3]
        bp_processed += args.window_size

        previous_update_time = progress_meter(previous_update_time, chrm, stop, bp_processed, total_bp_in_dataset)


if __name__ == '__main__':
    main()