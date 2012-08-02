#!/usr/bin/env python
# encoding: utf-8

"""
run_dadi.py

Created by Nick Crawford on 2012-07-05.
Copyright (c) 2012

The author may be contacted at ngcrawford@gmail.com


python VCF2Dadi.py \
-i test_data/butterfly.vcf.gz \
-o test_data/butterfly.dadi.input.txt \
-L 1:500-1000 \
-w 5000 \
-overlap 0 \
-p cydno:c511,c512,c513,c514,c515,c563,c614,c630,c639,c640 \
outgroups:h665,i02-210 \
melpo:m523,m524,m525,m589,m675,m676,m682,m683,m687,m689 
"""

import sys
import VCF
import dadi
import types
import argparse
import itertools
import multiprocessing
from copy import copy, deepcopy

def get_args():
    """Parse sys.argv"""
    parser = argparse.ArgumentParser()

    parser.add_argument('-i','--input', required=True,
                        help='Path to VCF file.')
    
    parser.add_argument('-o','--output',
                        help='Path to output csv file. If path is not set defaults to STDOUT.')

    parser.add_argument('-p','--populations', nargs='+',
                        help='Names of populations and samples. The format is: "PopName:sample1,sample2,sample3,etc..."')

    parser.add_argument('-L','--region', default=None, type=str,
                        help='chrm:start-stop')

    parser.add_argument('-w', '--window-size', type=int,
                        help='The size of the windows')

    parser.add_argument('-overlap', type=int, default=0,
                        help='The number of base pairs each window overlaps the previous')

    args = parser.parse_args()

    populations_dict  = {}
    for pop in args.populations:
        pop_name, sample_ids = pop.strip().split(":")
        sample_ids = sample_ids.split(",")
        populations_dict[pop_name] = sample_ids

    args.populations = populations_dict

    if args.region != None:
        if len(args.region.split(":")) == 2:
            chrm = [args.region.split(":")[0]]
            start_stop = [int(item) for item in args.region.split(":")[1].split("-")]
            args.region = chrm + start_stop

    else:
        args.region = [args.region]


    return args


def create_dadi_header(args):
    pop_ids = args.populations.keys()
    dadi_header = ['Outgroup','Ingroup','Allele1','Allele2','Chrm','Pos']
    dadi_header[3:3] = pop_ids
    dadi_header[-2:2] = pop_ids
    dadi_header = ' '.join(dadi_header)
    return dadi_header


def create_Fstats_header(pop_ids):

    pop_values = ['Tajimas_D','W_theta','pi','Seg_Sites']

    final_header = ['chrm', 'start','stop','Fst']
    for count, pop in enumerate(pop_ids):
        final_header += [pop + "." + i for i in pop_values]

    return final_header

def target(args):
    """Takes an object and a list of arguments the first
       of which is the method name to run.
    """
    object = args[0]
    method_name = args[1]
    return getattr(object, method_name)(*args[2:])

def make_vcf_slices(slices, vcf, args):
    for ccount, chrm in enumerate(slices.keys()):
        for scount, s in enumerate(slices[chrm]):
            yield tuple([vcf, vcf.slice_vcf.__name__, args.input, chrm] + list(s))

def process_vcf_slices(slices, vcf):        
    for s in slices:
        if s != []:
            yield tuple([vcf, vcf.count_alleles_in_vcf.__name__] + [(s)])


def calc_fstats_with_dadi():
    """Takes a VCF file and calculates """

    vcf = VCF.VCF()
    vcf.set_header(args.input)
    vcf.populations = args.populations
    slices = vcf.generate_slices(args)
    pool = multiprocessing.Pool(2)

    sliced_vcf = itertools.imap(target, make_vcf_slices(slices, vcf, args))

    chunksize = 1000000 * args.window_size # if 1 snps per kb then this works out to ~ 
    processed_vcf = pool.imap(target, process_vcf_slices(sliced_vcf, vcf))
    
    fout = open(args.output,'w')
    header = ' '.join(create_Fstats_header(vcf.populations)) + '\n'
    fout.write(header)

    for count, s in enumerate(processed_vcf):
        region = s[0]
        region = [region['CHROM'], region['POS']]

        dd = vcf.make_dadi_fs(s)
        
        if dd == None: continue # skip empty calls
        
        pop_ids = vcf.populations.keys()
        pop_ids.remove('outgroups')
        projection_size = 10
        pairwise_fs  = dadi.Spectrum.from_data_dict(dd, pop_ids, [projection_size]*2)

        # Create final line, add Fst info
        final_line = region
        final_line += [pairwise_fs.Fst()]

        # Add in population level stats
        for pop in pop_ids:
            fs = dadi.Spectrum.from_data_dict(dd, [pop], [projection_size])
            final_line += [fs.Tajima_D(), fs.Watterson_theta(), fs.pi(), fs.S()]

        # write output
        final_line = [str(i) for i in final_line]
        fout.write(' '.join(final_line) + "\n")

    fout.close()
        
if __name__ == '__main__':
    args = get_args()
    calc_fstats_with_dadi()



