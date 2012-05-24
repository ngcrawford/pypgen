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

def _pickle_method(method):
    func_name = method.im_func.__name__
    obj = method.im_self
    cls = method.im_class

    return _unpickle_method, (func_name, obj, cls)

def _unpickle_method(func_name, obj, cls):
    for cls in cls.mro():
        try:
            func = cls.__dict__[func_name]
        except KeyError:
            pass
        else:
            break
    return func.__get__(obj, cls)

class VCF(object):
    """docstring for VCF"""
    def __init__(self):
        super(VCF, self).__init__()
    
        self.header = None
        self.__header_dict__ = None

        self.sample_format = None
        self.sample_format_dict = None

        self.populations = None
        self.f_statistic  = None
        self.vcf_file = None
        self.window_size = None

    def parse_individual_snps(self, vcf, fout, stat_id):

        vcf = open(vcf,'rU')
        
        if fout == None:
            fout = sys.stdout
        else:
            fout = open(fout,'w')

        # SETUP NAMED TUPLE TO STORE INFO FROM A SINGLE BASE
        field_labels = []

        # SETUP COUNTERS
        current_base = None
        line_count = None

        # PARSE VCF FIlE
        snp_count = 0

        for line in vcf:

            # SKIP HEADER
            if line.startswith("#CHROM"):
                self.header = line.strip("#").strip().split()
                self.__header_dict__ = OrderedDict([(item,None) for item in self.header])
                line_count = 0

            # START PROCESSING ALIGNED BASES
            if line_count > 0:
                vcf_line_dict = self.parse_vcf_line(line)
  
                allele_counts = self.calc_allele_counts(vcf_line_dict)
                f_statistics = self.calc_fstats(allele_counts)
                info = self.parse_info(vcf_line_dict["INFO"])
                formated_data, header = self.__write_to_outfiles__(vcf_line_dict["CHROM"],\
                                    vcf_line_dict["POS"],info["DP"], stat_id, f_statistics)

                if line_count == 1:
                    fout.write(header+"\n")
                    fout.write(formated_data+"\n")
                else:
                    fout.write(formated_data+"\n")

            if line_count >= 0:
                line_count += 1

    def parse_info(self,info_field):
        
        info = []
        for item in info_field.split(','): # TO DO: comma should be ';'
            pair = item.split("=") 
            if len(pair) == 2:
                info.append(pair)

        info_dict = dict(info)
        return info_dict

    def parse_vcf_line(self, line):

        line_parts = line.strip().split()
        vcf_line_dict = self.__header_dict__

        for count, item in enumerate(vcf_line_dict):
            vcf_line_dict[item] = line_parts[count]

        if self.sample_format == None:
            self.sample_format = vcf_line_dict["FORMAT"].split(":")
            # self.sample_format_dict = dict([(item,None) for item in self.sample_format])      
     
        for count, item in enumerate(vcf_line_dict):
            
            if count >= 9:
                genotype = vcf_line_dict[item]
                
                if genotype == "./.":
                    vcf_line_dict[item] = None

                else:
                    genotype = dict(zip(self.sample_format,genotype.split(":")))
                    vcf_line_dict[item] = genotype

        return vcf_line_dict

    def calc_allele_counts(self, vcf_line_dict):

        allele_counts = self.populations.fromkeys(self.populations.keys(),None)

        for population in self.populations.keys():

            allele_format_dict = {0:0,1:0,2:0,3:0,4:0}   # create dict to prevent pointer issues
            allele_counts[population] = allele_format_dict

            for sample_id in self.populations[population]:
                

                if vcf_line_dict[sample_id] != None:
                    
                    genotype = vcf_line_dict[sample_id]
                    genotype = genotype["GT"].split("/")

                    if genotype == [".","."]: continue

                    genotype = [int(item) for item in genotype]
                    
                    for allele in genotype:
                        allele_counts[population][allele] += 1 
        
        return allele_counts

   
    def calc_fstats(self, allele_counts):

        #test_pair = {'So_Riviere_Goyaves': np.array([ 0.0, 1.0, 0.0, 0.0]), 'Plage_de_Viard': np.array([ 1.0, 0.0, 0.0, 0.0]),}
        
        # CALCULATE ALLELE FREQUENCIES
        allele_freqs_dict = self.populations.fromkeys(self.populations.keys(),None)
        for population in allele_counts.keys():
            counts =  allele_counts[population].values()
            freqs =  counts/np.sum(counts,dtype=float)
            allele_freqs_dict[population] = freqs

        allele_freqs = allele_freqs_dict.values()

        # CACULATE PAIRWISE F-STATISTICS
        pairwise_results = {}
        for population_pair in combinations(self.populations.keys(),2):
            
            pop1, pop2 =  population_pair
            Ns = [sum(allele_counts[pop].values()) for pop in [pop1, pop2]]

            if 0 in Ns: 
                values = [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]
                values_dict = dict(zip(['Hs_est', 'Ht_est', 'Gst_est', 'G_prime_st_est', 'G_double_prime_st_est', 'D_est'],values))
                pairwise_results[population_pair] = values_dict
                continue
            
            else:

                pop1 = allele_freqs_dict[pop1]
                pop2 = allele_freqs_dict[pop2]
                allele_freqs = [pop1, pop2]

                n = float(len(Ns))
                Ns_harm = fstats.harmonic_mean(Ns)
                Ns_harm_chao = fstats.harmonic_mean_chao(Ns)

                # CALCULATE Hs AND Ht
                n = 2
                Hs_prime_est_ = fstats.Hs_prime_est(allele_freqs,n)
                Ht_prime_est_ = fstats.Ht_prime_est(allele_freqs,n)
                Hs_est_ = fstats.Hs_est(Hs_prime_est_, Ns_harm)
                Ht_est_ = fstats.Ht_est(Ht_prime_est_,Hs_est_,Ns_harm,n)

                # CALCULATE F-STATISTICS
                Gst_est_ = fstats.Gst_est(Ht_est_, Hs_est_)
                G_prime_st_est_ = fstats.G_prime_st_est(Ht_est_, Hs_est_, Gst_est_, n)
                G_double_prime_st_est_ = fstats.G_double_prime_st_est(Ht_est_, Hs_est_, n)
                D_est_ = fstats.D_est(Ht_est_, Hs_est_, n)
                
                # PRINT OUTPUT
                values = [Hs_est_, Ht_est_, Gst_est_, G_prime_st_est_, G_double_prime_st_est_, D_est_]
                values_dict = dict(zip(['Hs_est', 'Ht_est', 'Gst_est', \
                                        'G_prime_st_est', 'G_double_prime_st_est', 'D_est'],values))

                pairwise_results[population_pair] = values_dict

        return pairwise_results

    def __write_to_outfiles__(self, chrm, pos, depth, stat_id, f_statistics):

        line = ",".join((chrm, pos, str(depth)))
        for pair in f_statistics.keys():
            if f_statistics[pair] == None:
                continue

            if f_statistics[pair][stat_id] == np.nan:
                value = 'NA'
            else:
                value = str(f_statistics[pair][stat_id])
            line += "," + value
        
        header = 'chrm,pos,total_depth,' + ','.join([left+"-"+right for left, right in f_statistics.keys()]) + "," + stat_id
        return (line, header)
     
    def set_header(self, vcf_path):
        vcf_file = open(vcf_path,'rU')
        for line in vcf_file:
            if line.startswith("#CHROM"):
                self.header = line.strip("#").strip().split()
                self.__header_dict__ = OrderedDict([(item,None) for item in self.header])

        vcf_file.close()
       

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





if __name__ == '__main__':

    args = get_args()
    copy_reg.pickle(types.MethodType, _pickle_method, _unpickle_method)
    vcf = VCF()
    vcf.vcf_file = args.input
    vcf.set_header(vcf.vcf_file)
    vcf.populations = args.populations
    vcf.f_statistic = args.f_statistic
    vcf.window_size = args.window_size

    if args.window_size != None:
        do_windowed_analysis(vcf)

    else:
        vcf.populations = args.populations
        vcf.parse_individual_snps(args.input, args.output, args.f_statistic)



