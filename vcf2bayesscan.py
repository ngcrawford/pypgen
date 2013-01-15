#!/usr/bin/env python
# encoding: utf-8

"""
[loci]=1000

[populations]=10

[pop]=1
  1  40  2    0  40
  2  40  2   34  6
  3  40  2   19  21
  4  40  2   31  9

"""

"""
python vcf2bayesscan.py \
-i test_data/butterfly.vcf.gz -o bayesscan.snps \
-p cydno:c511,c512,c513,c514,c515,c563,c614,c630,c639,c640 \
outgroups:h665,i02-210 \
melpo:m523,m524,m525,m589,m675,m676,m682,m683,m687,m689 \
pachi:p516,p517,p518,p519,p520,p591,p596,p690,p694,p696 

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
from collections import defaultdict


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

    sample_ids =  [s for p in populations.values() for s in p]

    output_population_info = defaultdict(list)

    snp_count = 0
    for count, line in enumerate(open_vcf(args)):
        
        if line.startswith('#') == True: continue
        snp_count += 1

        vcf_line = parse_vcf_line(line, empty_vcf_line)
        allele_counts = calc_allele_counts(populations, vcf_line)

        possible_alleles = set([a for d in allele_counts.keys() \
                                  for a in allele_counts[d].keys() \
                                  if allele_counts[d][a] > 0.0 ])

        
        for pop in populations.keys():
            genes = len([vcf_line[s] for s in populations[pop] if vcf_line[s] != None]) * 2
            info = '\t'.join(map(str, [snp_count, genes, len(possible_alleles)]))
            a_counts = "\t".join([str(int(allele_counts[pop][a])) for a in possible_alleles])
            line = "{}\t{}\n".format(info, a_counts)
            output_population_info[pop].append(line)


    fout = open('test.bayesscan.output','w')
    fout.write("[loci]={}\n\n".format(len(output_population_info[output_population_info.keys()[0]])))
    fout.write("[populations]={}\n\n".format(len(populations.keys())))

    for count, pop in enumerate(output_population_info.keys()):

        fout.write('[pop]={}\n'.format(pop))
        
        for line in output_population_info[pop]:
            fout.write(line)






if __name__ == '__main__':
    main()






