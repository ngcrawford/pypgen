#!/usr/bin/env python
# encoding: utf-8

import os
import re
import sys
import gzip
import pysam
import itertools
import numpy as np
from fstats import fstats
from collections import OrderedDict
from itertools import combinations, izip_longest

def set_header(vcf_path):
    """Open VCF file and read in header line as Ordered Dict"""

    vcf_file = gzip.open(vcf_path,'rb')
    header_dict = None
    for line in vcf_file:
        if line.startswith("#CHROM"):
            header = line.strip("#").strip().split()
            header_dict = OrderedDict([(item, None) for item in header])
            break

    vcf_file.close()
    return header_dict

def parse_populations_list(populations):
    populations_dict  = {}
    for pop in populations:
        pop_name, sample_ids = pop.strip().split(":")
        sample_ids = sample_ids.split(",")
        populations_dict[pop_name] = sample_ids

    return populations_dict

def pairwise(iterable):
    """Generates paris of slices from iterator
       s -> (s0,s1), (s1,s2), (s2, s3), ..."""

    a, b = itertools.tee(iterable)
    next(b, None)
    return itertools.izip(a, b)

def get_slice_indicies(vcf_bgzipped_file, window_size=500):
    """Get slice information from VCF file that is tabix indexed file (bgzipped). """

    # READ IN FILE HEADER
    tbx = pysam.Tabixfile(vcf_bgzipped_file) # TODO: create try statement to test that file is actually a VCF
    chrms = tbx.contigs

    # PARSE LENGTH INFO FROM HEADER
    chrm_lengths = []
    for line in tbx.header:
        
        if line.startswith("##contig="):
            
            chrm_name = re.findall(r'ID=.*,',line)
            chrm_name = chrm_name[0].strip('ID=').strip(',')
            
            chrm_length = re.findall(r'length=.*>',line)
            chrm_length = int(chrm_length[0].strip('length=').strip('>'))
            
            chrm_lengths.append((chrm_name, 1, chrm_length))

    chrm_lengths = tuple(chrm_lengths)
    tbx.close()

    # GENERATE SLICES
    # current does not make overlapping slices. 
    for chrm, start, stop in chrm_lengths:

        slice_indicies = itertools.islice(xrange(start, stop + 1), 0, stop + 1, window_size)
        
        for count, s in enumerate(pairwise(slice_indicies)):
            yield (chrm, s[0], s[1] -1) # subtract one to prevent 1 bp overlap

def slice_vcf(vcf_bgzipped_file, chrm, start, stop):
    
    tbx = pysam.Tabixfile(vcf_bgzipped_file)
    return tuple(row for row in tbx.fetch(chrm, start, stop))


def parse_info_field(info_field):
    
    info_dict = {}
    for item in info_field.split(';'):
        pair = item.split("=") 
        if len(pair) == 2:
            info_dict[pair[0]] = pair[1] # this could be improved on
    return info_dict

def parse_vcf_line(pos, header):
    """Read in VCF line and convert it to an OrderedDict"""

    pos_parts = pos.strip().split()
    vcf_line_dict = header

    for count, item in enumerate(vcf_line_dict):
        vcf_line_dict[item] = pos_parts[count]

    sample_format = vcf_line_dict["FORMAT"].split(":")
 
    for count, item in enumerate(vcf_line_dict):
        
        if count >= 9:
            genotype = vcf_line_dict[item]
            
            if genotype == "./." or genotype == ".":      # "'./.'' for dip, '.' for haploid
                vcf_line_dict[item] = None

            else:
                genotype = dict(zip(sample_format,genotype.split(":")))

                # CONVERT STRINGS TO APPOPRIATE TYPES (INTS, FLOATS, ETC.)
                genotype['GQ'] = float(genotype['GQ'])
                genotype['DP'] = int(genotype['DP'])
                genotype['AD'] = tuple(int(ad) for ad in genotype['AD'].split(","))
                genotype['PL'] = tuple(int(ad) for ad in genotype['PL'].split(","))

                vcf_line_dict[item] = genotype
    
    vcf_line_dict['POS'] = int(vcf_line_dict['POS'])
    vcf_line_dict['QUAL'] = float(vcf_line_dict['QUAL'])

    vcf_line_dict['INFO'] = parse_info_field(vcf_line_dict['INFO'])

    return vcf_line_dict




def filter_on_population_sizes(vcfline, populations, min_samples=0):

    sample_counts = []

    for pop in populations.keys():
        sample_count = 0
        for sample in populations[pop]:
            if vcfline[sample] != None:
                sample_count += 1

        sample_counts.append(sample_count)

    sample_counts = [item for item in sample_counts if item >= min_samples]

    if len(sample_counts) >= 2: 
        return True
    else: 
        return False


def format_output(chrm, start, stop, depth, stat_id, multilocus_f_statistics):

    line = ",".join((chrm, str(start), str(stop), str(depth)))
    for pair in multilocus_f_statistics.keys():
        if multilocus_f_statistics[pair] == None:
            continue

        if multilocus_f_statistics[pair][stat_id] == np.nan:
            value = 'NA'
        
        else:
            value = str(multilocus_f_statistics[pair][stat_id])
        
        line += "," + value
    
    line += "\n"
    return (line)

def calc_allele_counts(populations, vcf_line_dict):

    allele_counts = populations.fromkeys(populations.keys(), None)

    for population in populations.keys():

        allele_format_dict = {0:0.0,1:0.0,2:0.0,3:0.0,4:0.0}   # create dict to prevent pointer issues
        allele_counts[population] = allele_format_dict

        for sample_id in populations[population]:

            if vcf_line_dict[sample_id] != None:
                
                genotype = vcf_line_dict[sample_id]
                genotype = genotype["GT"].split("/")

                if genotype == [".","."]: continue

                genotype = [int(item) for item in genotype]
                
                for allele in genotype:
                    allele_counts[population][allele] += 1.0 
    
    return allele_counts

def calc_fstats(populations, allele_counts):


    #test_pair = {'So_Riviere_Goyaves': np.array([ 0.0, 1.0, 0.0, 0.0]), 'Plage_de_Viard': np.array([ 1.0, 0.0, 0.0, 0.0]),}
    
    # CALCULATE ALLELE FREQUENCIES
    allele_freqs_dict = populations.fromkeys(populations.keys(), None)
    for population in allele_counts.keys():
        counts =  allele_counts[population].values()
        
        if sum(counts) == 0.0:
            freqs = [0.0] * 4

        else:
            freqs =  [count/sum(counts) for count in counts]

        allele_freqs_dict[population] = freqs

    allele_freqs = allele_freqs_dict.values()


    # CACULATE PAIRWISE F-STATISTICS
    pairwise_results = {}
    for population_pair in combinations(populations.keys(),2):
        
        pop1, pop2 =  population_pair
        Ns = [sum(allele_counts[pop].values()) for pop in [pop1, pop2]]

        if 0 in Ns: 
            values = [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]
            values_dict = dict(zip(['Hs_est', 'Ht_est', 'Gst_est', 
                                    'G_prime_st_est', 'G_double_prime_st_est', 
                                    'D_est'],
                                    values))
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
            Hs_prime_est_ = fstats.Hs_prime_est(allele_freqs, n)
            Ht_prime_est_ = fstats.Ht_prime_est(allele_freqs, n)
            Hs_est_ = fstats.Hs_est(Hs_prime_est_, Ns_harm)
            Ht_est_ = fstats.Ht_est(Ht_prime_est_, Hs_est_, Ns_harm, n)

            # CALCULATE F-STATISTICS
            Gst_est_ = fstats.Gst_est(Ht_est_, Hs_est_)
            G_prime_st_est_ = fstats.G_prime_st_est(Ht_est_, Hs_est_, Gst_est_, n)
            G_double_prime_st_est_ = fstats.G_double_prime_st_est(Ht_est_, Hs_est_, n)
            D_est_ = fstats.D_est(Ht_est_, Hs_est_, n)
            
            # PRINT OUTPUT
            values = [Hs_est_, Ht_est_, Gst_est_, G_prime_st_est_, 
                      G_double_prime_st_est_, D_est_]
            values_dict = dict(zip(['Hs_est', 'Ht_est', 'Gst_est', \
                                    'G_prime_st_est', 'G_double_prime_st_est', 
                                    'D_est'],
                                    values))

            pairwise_results[population_pair] = values_dict

    return pairwise_results

def output_header():
    return 'Chrm,Start,Stop,DP_mean,DP_std'

def calc_slice_stats(data):
    """Main function for caculating statistics.

       Make it easy to add more statitics.
    """

    tabix_slice, chrm, start, stop, populations, header = data

    Hs_est_list = []
    Ht_est_list = []

    chrm = None
    pos = None
    total_depth = []

    output_line = [chrm, start, stop]

    for count, line in enumerate(tabix_slice):

        vcf_line_dict = parse_vcf_line(line, header)

        # CREATE FILTERS HERE:
        if filter_on_population_sizes(vcf_line_dict, populations, min_samples=5) == False:
            continue

        total_depth.append(int(vcf_line_dict['INFO']["DP"]))

    total_depth = np.array(total_depth)
    
    # UPDATE OUTPUT LINE
    output_line += [total_depth.mean(), total_depth.std()] 

    allele_counts = calc_allele_counts(populations, vcf_line_dict)
    f_statistics = calc_fstats(populations, allele_counts)
    print f_statistics[('melpo', 'pachi')]



    # return format_output(chrm, start, stop, mean_depth, f_statistic, multilocus_f_statistics)



