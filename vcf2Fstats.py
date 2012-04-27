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
import argparse
import numpy as np
from fstats import fstats
from itertools import combinations
from collections import OrderedDict

def get_args():
    """Parse sys.argv"""
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input', required=True,
                        help='Path to VCF file.')

    parser.add_argument('-o','--output', required=True,
                        help='Path to output csv file.')

    args = parser.parse_args()
    return args


class VCF(object):
    """docstring for VCF"""
    def __init__(self):
        super(VCF, self).__init__()
    
        self.header = None
        self.__header_dict__ = None

        self.sample_format = None
        self.sample_format_dict = None

        self.populations = None

    def set_populations(self, pop_dict):
        self.populations = pop_dict

    def parse_file_path(self, vcf, fout):

        vcf = open(vcf,'rU')
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
                fstats = self.fstats(allele_counts)
                info = self.parse_info(vcf_line_dict["INFO"])
                formated_data, header = self.__write_to_outfiles__(vcf_line_dict["CHROM"], vcf_line_dict["POS"],info["DP"],fstats)

                if line_count == 1:
                    fout.write(header+"\n")
                    fout.write(formated_data+"\n")
                else:
                    fout.write(formated_data+"\n")

            if line_count >= 0:
                line_count += 1

            # if line_count > 5: break


    def parse_info(self,info_field):
        
        info = []
        for item in info_field.split(';'):
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

            allele_format_dict = {0:0,1:0,2:0,3:0}   # create dict to prevent pointer issues
            allele_counts[population] = allele_format_dict

            for sample_id in self.populations[population]:
                
                if vcf_line_dict[sample_id] != None:
                    
                    genotype = vcf_line_dict[sample_id]
                    genotype = genotype["GT"].split("/")
                    genotype = [int(item) for item in genotype]
                    
                    for allele in genotype:
                        allele_counts[population][allele] += 1 
        
        return allele_counts

   
    def fstats(self, allele_counts):

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
        for population_pair in combinations(populations,2):
            
            pop1, pop2 =  population_pair
            Ns = [ sum(allele_counts[pop].values()) for pop in [pop1, pop2]]

            if 0 in Ns: 
                values = [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]
                pairwise_results[population_pair] = values
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
                values = [ Hs_est_, Ht_est_, Gst_est_, G_prime_st_est_, G_double_prime_st_est_, D_est_,]
                pairwise_results[population_pair] = values

        return pairwise_results

    def __write_to_outfiles__(self, chrm, pos, depth,fstats):

        line = ",".join((chrm, pos, str(depth)))
        for pair in fstats.keys():
            if fstats[pair][0] == np.nan:
                value = 'NA'
            else:
                value = str(fstats[pair][0])
            line += "," + value
        
        header = 'chrm,pos,total_depth,' + ','.join([left+"-"+right for left, right in fstats.keys()])
        return (line, header)

if __name__ == '__main__':

    populations = {"Capasterre": ['CEJ035', 'CEJ036', 'CEJ037', 'CEJ039',  'CEJ040', 'CEJ041',],
                "Marina": ['CEJ084', 'CEJ085', 'CEJ088', 'CEJ106', 'CEJ108'],
                "No_of_Goyave": ['CEJ111','CEJ112', 'CEJ119','CEJ120', 'CEJ121', 'CEJ122',],
                "Plage_au_Petit_Bourg": ['CEJ092', 'CEJ094', 'CEJ114', 'CEJ116', 'CEJ117'], 
                "Plage_de_Viard": ['CJS1974','CJS1976'],
                "So_Riviere_Goyaves": ['CJS2068','CJS2069','CJS2073', 'CJS2074'],
                "St_Marie": ['CEJ029', 'CEJ030', 'CEJ031','CEJ032','CEJ033', 'CEJ034']}

    args = get_args()
    vcf = VCF()
    vcf.set_populations(populations)
    vcf.parse_file_path(args.input, args.output)


