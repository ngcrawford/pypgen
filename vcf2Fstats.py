#!/usr/bin/env python
# encoding: utf-8

"""
run_dadi.py

Created by Nick Crawford on 2012-07-05.
Copyright (c) 2012

The author may be contacted at ngcrawford@gmail.com


python vcf2Fstats.py \
-i test_data/butterfly.vcf.gz \
-o test_data/butterfly.dadi.input.txt \
-L Chr01:0-5000 \
-w 500 \
-overlap 0 \
-n 3 \
-p cydno:c511,c512,c513,c514,c515,c563,c614,c630,c639,c640 \
outgroups:h665,i02-210 \
melpo:m523,m524,m525,m589,m675,m676,m682,m683,m687,m689 \
pachi:p516,p517,p518,p519,p520,p591,p596,p690,p694,p696

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

    parser.add_argument('-n', '--processors', type=int, default=0,
                        help='The number of processors to use.')

    parser.add_argument('--projection-size', type=int, default=10,
                        help='Number of alleles dadi project the data down to.')

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

def slices_2_calls(slices, vcf, args):
    for ccount, chrm in enumerate(slices.keys()):
        
        if slices[chrm] == None: 
            print 'skipping' ,chrm
            continue
        
        for scount, s in enumerate(slices[chrm]):
            yield tuple([vcf, vcf.slice_2_allele_counts.__name__, args.input, chrm] + list(s))



def calc_fstats_with_dadi(args):
    """Takes a VCF file and calculates """

    vcf = VCF.VCF()
    print 'Setting header...'
    vcf.set_header(args.input)
    vcf.populations = args.populations
    pool = multiprocessing.Pool(args.processors)
    
    print 'Generating slices...'
    slices = vcf.generate_slices(args)
        
    # Creating output
    fout = open(args.output,'w')
    header = ' '.join(create_Fstats_header(vcf.populations)) + '\n'
    fout.write(header)

    z = slices_2_calls(slices, vcf, args)

    current_chrm = None
    for count, s in enumerate(pool.imap(target, slices_2_calls(slices, vcf, args), chunksize=1)):

        s = [i for i in s if i != None]
        region = [s[0]['CHROM'], s[0]['POS']]

        dd = vcf.make_dadi_fs(s)
        
        if current_chrm == None or current_chrm != s[0]['CHROM']:
            print 'Processing:', s[0]['CHROM']
            current_chrm = s[0]['CHROM']

        if dd == None: continue # skip empty calls
        
        pop_ids = vcf.populations.keys()
        pop_ids.remove('outgroups')
        projection_size = 10
        pairwise_fs  = dadi.Spectrum.from_data_dict(dd, pop_ids, [projection_size]*len(pop_ids))
        try:
            pairwise_fs  = dadi.Spectrum.from_data_dict(dd, pop_ids, [projection_size]*len(pop_ids))
        except:
            continue

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
    calc_fstats_with_dadi(args)



