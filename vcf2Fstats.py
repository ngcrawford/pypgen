#!/usr/bin/env python
# encoding: utf-8

"""
vcf2oneliners.py

Created by Nick Crawford on 2012-24-04.
Copyright (c) 2012

The author may be contacted at ngcrawford@gmail.com
"""

import os
import sys
import math
import argparse
import copy_reg
import types
from VCF import VCF
import multiprocessing
import numpy as np
from fstats import fstats
from itertools import combinations, izip_longest
from collections import OrderedDict

def get_args():
    """Parse sys.argv"""
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input', required=True,
                        help='Path to VCF file.')
    parser.add_argument('-o','--output',
                        help='Path to output csv file. If path is not set defaults to STDOUT.')
    parser.add_argument('-c','--cores', required=True,
                        help='Number of cores to use.')
    parser.add_argument('-p','--populations', nargs='+',
                        help='Names of populations and samples. The format is: "PopName:sample1,sample2,sample3,etc..."')

    parser.add_argument('-w','--window-size',default=None, type=int,
                        help='Size of the window in which to calculate pairwise F-staticstics')
    parser.add_argument("-f",'--f-statistic',required=True, 
                        choices=['Gst_est', 'G_prime_st_est', 'G_double_prime_st_est', 'D_est'])

    args = parser.parse_args()

    populations_dict  = {}
    for pop in args.populations:
        pop_name, sample_ids = pop.strip().split(":")
        sample_ids = sample_ids.split(",")
        populations_dict[pop_name] = sample_ids

    args.populations = populations_dict

    return args

def grouper(n, iterable, padvalue=None):
    """grouper(3, 'abcdefg', 'x') -->
    ('a','b','c'), ('d','e','f'), ('g','x','x')"""

    return izip_longest(*[iter(iterable)]*n, fillvalue=padvalue)


def slidingWindow(vcf, window_size=1000):
    """Generator function that yields non overlapping chunks of VCF lines."""

    # TO DO: add overlapping increment
    chrm_id = None
    start = 0
    stop = window_size
    chunk = []

    for count, line in enumerate(vcf):
        
        # SKIP HEADER
        if line.startswith("#"):
            continue

        # GET CHRM AND POS
        line_parts = line.strip().split()
        chrm, pos = line_parts[:2]
        pos = int(pos)

        # UPDATE CHRM ID FOR INITAL SITE
        if chrm_id == None:
            chrm_id = chrm

        # REZERO AT NEW CHRM
        #   AND YIELD CHUNK
        if chrm_id != chrm:
            yield chunk
            stop = window_size 
            chrm_id = chrm
            chunk = [line]
            continue

        # UPDATE CHUNK
        if pos < stop:
            chunk.append(line)

        # YIELD CHUNK IF CRURRENT 
        #   POS EXCEEDS STOP
        if pos >= stop:
            yield chunk
            stop += window_size
            chunk = [line]


def process_window(window):
    vcf = VCF()

    Hs_est_list = []
    Ht_est_list = []

    chrm = None
    pos = None
    depth = None

    for count, line in enumerate(window):
        vcf_line_dict = vcf.parse_vcf_line(line)
        print vcf_line_dict
        if count == 0:
            chrm = vcf_line_dict["CHROM"]
            pos = vcf_line_dict['POS']
            info = vcf.parse_info(vcf_line_dict["INFO"])
            depth = info["DP"]

        allele_counts = vcf.calc_allele_counts(vcf_line_dict)
        f_statistics = vcf.calc_fstats(allele_counts)

        Hs_est = f_statistics[('MAR', 'CAP')]['Hs_est']
        Hs_est_list.append(Hs_est)

        Ht_est = f_statistics[('MAR', 'CAP')]['Ht_est']
        Ht_est_list.append(Ht_est)

    print Ht_est_list

    # # Remove uninformative SNPs (nans)
    # Hs_est_list = [item for item in Hs_est_list if not math.isnan(item)]  
    # Ht_est_list = [item for item in Ht_est_list if not math.isnan(item)]

    # multilocus_f_statistics = {('MAR', 'CAP'): None}
    # if len(Hs_est_list) != 0 or len(Ht_est_list)!= 0:

    #     n = 2 # fix this
    #     Gst_est = fstats.multilocus_Gst_est(Ht_est_list, Hs_est_list)
    #     G_prime_st_est = fstats.multilocus_G_prime_st_est(Ht_est_list, Hs_est_list, n)
    #     G_double_prime_st_est = fstats.multilocus_G_double_prime_st_est(Ht_est_list, Hs_est_list, n)
    #     D_est = fstats.multilocus_D_est(Ht_est_list, Hs_est_list, n)

    #     values_dict = dict(zip(['Gst_est', 'G_prime_st_est', 'G_double_prime_st_est', 'D_est'],\
    #                   [Gst_est, G_prime_st_est, G_double_prime_st_est, D_est]))

    #     multilocus_f_statistics[('MAR', 'CAP')] = values_dict

    # return vcf.__write_to_outfiles__( chrm, pos, depth, 'Gst_est', multilocus_f_statistics)


def do_windowed_analysis(vcf):

    fout = open(args.output, 'w')
    vcf_file = open(vcf.vcf_file,'rU')

    p = multiprocessing.Pool(4)

    for chunk in grouper(10, slidingWindow(vcf_file,vcf.window_size)):
        results = p.map(process_window, chunk)
        for r in results:
            print r     

        # for count, window in enumerate(vcf.slidingWindow(vcf_file,args.window_size)):
        #     formated_data, header = vcf.process_window(window, args)
            
        #     if line_count == 0:
        #         fout.write(header+"\n")
        #         fout.write(formated_data+"\n")
        #     else:
        #         fout.write(formated_data+"\n")
        # fout.close()

        return self.__write_to_outfiles__( chrm, pos, depth, args.f_statistic, multilocus_f_statistics)

if __name__ == '__main__':

    args = get_args()
    copy_reg.pickle(types.MethodType, _pickle_method, _unpickle_method)
    vcf = VCF()
    vcf.vcf_file = args.input
    vcf.set_header(vcf.vcf_file)
    vcf.populations = args.populations
    vcf.f_statistic = args.f_statistic
    vcf.window_size = args.window_size





