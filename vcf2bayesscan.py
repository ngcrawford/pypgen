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

import os
import sys
import argparse
from VCF import *
import multiprocessing

"""
python vcf2bayesscan.py \
-i test_data/butterfly.vcf.gz -o bayesscan.snps \
-p cydno:c511,c512,c513,c514,c515,c563,c614,c630,c639,c640 \
outgroups:h665,i02-210 melpo:m523,m524,m525,m589,m675,m676,m682,m683,m687,m689 \
pachi:p516,p517,p518,p519,p520,p591,p596,p690,p694,p696 

"""



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

    parser.add_argument("-m",'--min-samples',
                        type=int,
                        default=5,
                        help="Minimum number of samples per population.")

    return parser.parse_args()

def call_genotype(gt_string):
    final_call = None
    call = gt_string.split(":")

def main():
    # get args. 
    args = get_args()
    populations = parse_populations_list(args.populations)
    print populations

    header = set_header(args.input)

    fin = gzip.open(args.input,'rb')

    outfiles = dict([(pop, open("{}_snps.txt".format(pop),'w')) for pop in populations.keys()])
    
    for count, line in enumerate(fin):
        if line.startswith('#') == True: continue
        if count > 10: break
        vcf_line_dict = parse_vcf_line(line, header)
        
        for pop in populations.keys():

            sample_count = 0
            calls = {'0':0, '1':0, '2':0, '3':0,}
            for sample in populations[pop]:
                if vcf_line_dict[sample] == None: continue

                gt =  vcf_line_dict[sample]["GT"]
                gt = gt.split('/')
                calls[gt[0]] += 1
                calls[gt[1]] += 1

                sample_count += 1

            print pop, count, sample_count, calls



if __name__ == '__main__':
    main()






