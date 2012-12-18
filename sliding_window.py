#!/usr/bin/env python
# encoding: utf-8

import os
import sys
import argparse
from VCF import *
import multiprocessing

class FooAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        print '%r %r %r' % (namespace, values, option_string)
        setattr(namespace, self.dest, values)

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
                        type=str,
                        help="Define a chromosomal region. \
                              A region can be presented, for example, in the following \
                              format: ‘chr2’ (the whole chr2), ‘chr2:1000000’ (region \
                              starting from 1,000,000bp) or ‘chr2:1,000,000-2,000,000’ \
                              (region between 1,000,000 and 2,000,000bp including the end \
                              points). The coordinate is 1-based.' [Same format as \
                              SAMTOOLs/GATK, example text cribbed from SAMTOOLs]")
    
    parser.add_argument('-p','--populations', 
                        nargs='+',
                        help='Names of populations and samples. \
                              The format is: "PopName:sample1,sample2,sample3,etc..."')

    parser.add_argument('-w','--window-size',
                        default=None, 
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
               chrm, start, stop, populations, header, args.min_samples]

        # if count > 1: break

def main():
    # get args. 
    args = get_args()
    print args

    # TODO: 
    # test that pysam is installed.
    # bgzip check. MDSum?

    # Get slice indicies
    # 1. read file and get chrm sizes
    # 2. process chrm sizes and return as
    #    slices and a zipped list (chrm, (start, stop))
    slice_indicies = get_slice_indicies(args.input, args.regions, args.window_size)
    
    # Get information about samples from the header.
    # this becomes the precursor to the VCF row
    header = set_header(args.input)
    populations = parse_populations_list(args.populations)

    p = multiprocessing.Pool(4)
    order = []

    for count, result in enumerate(p.map(calc_slice_stats, generate_fstats_from_vcf_slices(slice_indicies, populations, header, args))):
        if result == None: continue
        chrm_start_stop, result = result
        stats, order = multilocus_f_statistics_2_sorted_list(result, order=order)

        if count == 0:
            print ','.join(['chrm','start','stop','snp_count', 'total_depth_mean', 'total_depth_stdev'] +  map(str, order))
            print ','.join(map(str, chrm_start_stop) + map(str, stats))
        else:
            print ','.join(map(str, chrm_start_stop) + map(str, stats))


if __name__ == '__main__':
    main()