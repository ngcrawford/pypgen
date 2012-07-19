import os
import sys
import math
import types
import argparse
import copy_reg
import numpy as np
import multiprocessing
from fstats import fstats
from collections import OrderedDict
from itertools import combinations, izip_longest



def get_args():
    """Parse sys.argv"""
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input', required=True,
                        help='Path to VCF file.')
    
    parser.add_argument('-o','--output',
                        help='Path to output csv file. If path is not set defaults to STDOUT.')
    
    parser.add_argument('-c','--cores', required=True, type=int,
                        help='Number of cores to use.')
    
    parser.add_argument('-p','--populations', nargs='+',
                        help='Names of populations and samples. The format is: "PopName:sample1,sample2,sample3,etc..."')

    parser.add_argument('-w','--window-size',default=None, type=int,
                        help='Size of the window in which to calculate pairwise F-staticstics')
    
    parser.add_argument("-f",'--f-statistic',required=True, 
                        choices=['Gst_est', 'G_prime_st_est', 'G_double_prime_st_est', 'D_est'])

    parser.add_argument("-m",'--min-samples',type=int,
                        help="Minimum number of samples per population.")

    args = parser.parse_args()

    populations_dict  = {}
    for pop in args.populations:
        pop_name, sample_ids = pop.strip().split(":")
        sample_ids = sample_ids.split(",")
        populations_dict[pop_name] = sample_ids

    args.populations = populations_dict

    return args


def parse_info(info_field):
    
    info = []
    for item in info_field.split(','): # TO DO: comma should be ';'
        pair = item.split("=") 
        if len(pair) == 2:
            info.append(pair)

    info_dict = dict(info)
    return info_dict

def parse_vcf_line(header_dict, line):

    line_parts = line.strip().split()
    vcf_line_dict = header_dict

    for count, item in enumerate(vcf_line_dict):
        vcf_line_dict[item] = line_parts[count]

    sample_format = vcf_line_dict["FORMAT"].split(":")
 
    for count, item in enumerate(vcf_line_dict):
        
        if count >= 9:
            genotype = vcf_line_dict[item]
            
            if genotype == "./.":
                vcf_line_dict[item] = None

            else:
                genotype = dict(zip(sample_format,genotype.split(":")))
                vcf_line_dict[item] = genotype

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
    if len(sample_counts) >= 2: return True
    else: return False

def calc_allele_counts(populations, vcf_line_dict):

    allele_counts = populations.fromkeys(populations.keys(),None)

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
            yield (chunk, stop-window_size, stop)
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
            yield (chunk, stop-window_size, stop)
            stop = pos - (pos % window_size)
            chunk = [line]

def slidingWindow_with_overlap(vcf, window_size=1000, overlap=500):
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
            yield (chunk, stop-window_size, stop)
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
            yield (chunk, stop-window_size, stop)
            stop = pos - (pos % window_size)
            chunk = [line]


 
def set_header(vcf_path):
    vcf_file = open(vcf_path,'rU')
    header_dict = None
    for line in vcf_file:
        if line.startswith("#CHROM"):
            header = line.strip("#").strip().split()
            header_dict = OrderedDict([(item, None) for item in header])
            break

    vcf_file.close()
    return header_dict

    vcf_file.close()

def grouper(n, iterable, padvalue=None):
    """grouper(3, 'abcdefg', 'x') -->
    ('a','b','c'), ('d','e','f'), ('g','x','x')"""

    return izip_longest(*[iter(iterable)]*n, fillvalue=padvalue)

def parse_populations_string(populations):
    populations_dict  = {}
    for pop in populations:
        pop_name, sample_ids = pop.strip().split(":")
        sample_ids = sample_ids.split(",")
        populations_dict[pop_name] = sample_ids

    return populations_dict

def process_window(data):

    window_parts, header_dict, populations, f_statistic = data

    if window_parts == None:
        return None
    else:
        window, start, stop = window_parts

    Hs_est_list = []
    Ht_est_list = []

    chrm = None
    pos = None
    depth = []

    for count, line in enumerate(window):
        vcf_line_dict = parse_vcf_line(header_dict, line)
        
        if count == 0:
            chrm = vcf_line_dict["CHROM"]
            pos = vcf_line_dict['POS']
            info = parse_info(vcf_line_dict["INFO"])
        

        if filter_on_population_sizes(vcf_line_dict, populations, min_samples=5) == False:
            continue

        depth.append(int(info["DP"]))

        allele_counts = calc_allele_counts(populations, vcf_line_dict)
        f_statistics = calc_fstats(populations, allele_counts)
        Hs_est = f_statistics[('MAR', 'CAP')]['Hs_est']
        Hs_est_list.append(Hs_est)

        Ht_est = f_statistics[('MAR', 'CAP')]['Ht_est']
        Ht_est_list.append(Ht_est)


    # Remove uninformative SNPs (nans)
    pairs = zip(Hs_est_list, Ht_est_list)
    pairs = [pair for pair in pairs if np.nan not in pair]

    
    multilocus_f_statistics = {('MAR', 'CAP'): None}
    if len(pairs) != 0:

        n = 2 # fix this
        Gst_est = fstats.multilocus_Gst_est(Ht_est_list, Hs_est_list)
        G_prime_st_est = fstats.multilocus_G_prime_st_est(Ht_est_list, Hs_est_list, n)
        G_double_prime_st_est = fstats.multilocus_G_double_prime_st_est(Ht_est_list, Hs_est_list, n)
        D_est = fstats.multilocus_D_est(Ht_est_list, Hs_est_list, n)

        values_dict = dict(zip(['Gst_est', 'G_prime_st_est', 'G_double_prime_st_est', 'D_est'],\
                      [Gst_est, G_prime_st_est, G_double_prime_st_est, D_est]))

        multilocus_f_statistics[('MAR', 'CAP')] = values_dict

    if len(depth) != 0:
         mean_depth = sum(depth)/float(len(depth))
    
    else:
        mean_depth = 0.0
        multilocus_f_statistics[('MAR', 'CAP')] = \
            dict(zip(['Gst_est', 'G_prime_st_est', 'G_double_prime_st_est', 'D_est'], ['nan']*4))

    return format_output(chrm, start, stop, mean_depth, f_statistic, multilocus_f_statistics)


def do_windowed_analysis(vcf_path, output, window_size, populations, f_statistic, cores):

    # OPEN INPUT AND OUTPUT FILES
    fout = open(output, 'w')
    vcf_file = open(vcf_path,'rU')
    
    # WRITE HEADER
    header_dict = set_header(vcf_path)
    left, right =populations.keys()
    fout.write('chrm,start,stop,mean_depth,' + ','.join([left+"-"+right]) + "-" + f_statistic+'\n')

    p = multiprocessing.Pool(cores)
    for chunk in grouper(5000, slidingWindow(vcf_file, window_size)):

        chunk = [(window, header_dict, populations, f_statistic) for window in chunk]
        results = p.map(process_window, chunk)

        for line in results:
            if line != None:
                fout.write(line)   

    fout.close()
    vcf_file.close()


def main():
    args = get_args()
    do_windowed_analysis(args.input, args.output, args.window_size, args.populations, args.f_statistic, args.cores)

if __name__ == '__main__':
    main()
