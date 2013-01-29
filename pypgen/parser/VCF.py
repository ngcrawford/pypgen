#!/usr/bin/env python
# encoding: utf-8

import os
import re
import sys
import gzip
import pysam
import argparse
import itertools
import numpy as np
from pypgen.misc import helpers
from pypgen.fstats import fstats
from collections import OrderedDict, defaultdict
from itertools import combinations, izip_longest


def default_args():
    """Parse sys.argv"""
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', '--input',
                        required=True,
                        type=str,
                        help='Path to VCF file.')

    parser.add_argument('-o', '--output',
                        nargs='?',
                        type=argparse.FileType('w'),
                        default=helpers.Unbuffered(sys.stdout),  # forces consistent writing to STOUT
                        help='Path to output csv file. \
                              If path is not set, defaults to STDOUT.')

    parser.add_argument('-c', '--cores',
                        required=True,
                        type=int,
                        help='Number of cores to use.')

    parser.add_argument('-r', '--regions',
                        required=False,
                        # action=FooAction, # TODO: fix this!
                        nargs='+',
                        help="Define a chromosomal region. \
                              A region can be presented, for example, in the following \
                              format: ‘chr2’ (the whole chr2), ‘chr2:1000000’ (region \
                              starting from 1,000,000bp) or ‘chr2:1,000,000-2,000,000’ \
                              (region between 1,000,000 and 2,000,000bp including the end \
                              points). The coordinate is 1-based.' Multiple regions can \
                              be submitted seperated by spaces. [NOte: this is the same \
                              format as SAMTOOLs/GATK, example text largely cribbed from \
                              SAMTOOLs]")

    parser.add_argument('-f', '--filter-string',
                        required=False,
                        default="FILTER == PASS")

    parser.add_argument('--regions-to-skip',
                        default=[],
                        required=False,
                        nargs='+',
                        help='Define a chromosomal region(s) to skip.')

    parser.add_argument('-p', '--populations',
                        nargs='+',
                        help='Names of populations and samples. \
                              The format is: "PopName:sample1,sample2,.. \
                              PopName2:sample3,sample4,..." \
                              Whitespace is used to delimit populations. Note: the \
                              population name uname "Outgroup" is reserved for \
                              samples that that are used to polarize genotype calls.')

    parser.add_argument('-w', '--window-size',
                        default=5000,
                        type=int,
                        help='Size of the window in which to \
                              calculate pairwise F-staticstics')

    parser.add_argument("-m", '--min-samples',
                        type=int,
                        default=5,
                        help="Minimum number of samples per population.")

    return parser


def process_snp_call(snp_call, ref, alt, IUPAC_ambiguities=False):
    """Process VCF genotype fields.
        The current version is very basic and
        doesn't directly take into account the
        quality of the call or call hets with 
        IUPAC ambiguity codes."""

    # IUPAC ambiguity codes
    IUPAC_dict = {('A','C'):'M',
                  ('A','G'):'R',
                  ('A','T'):'W',
                  ('C','G'):'S',
                  ('C','T'):'Y',
                  ('G','T'):'K',
                  ('A','C','G'):'V',
                  ('A','C','T'):'H',
                  ('A','G','T'):'D',
                  ('C','G','T'):'B'}

    #called_base = ""
    snp_call = snp_call.split(":")

    # process blanks
    if snp_call[0] == "./.":
        called_base = "-"

    else:
        allele1, allele2 = snp_call[0].split("/")

        # process "0/0"
        if allele1 == '0' and allele2 == '0':
            called_base = ref

        if allele1 == '1' and allele2 == '1':
            called_base = alt

        # process "0/N"
        if allele1 == '0' and allele2 != '0':

            if IUPAC_ambiguities == False:
                called_base = 'N'

            else:
                call = [ref] + [alt.split(',')[int(allele2) - 1]]
                call.sort()
                call = tuple(call)
                called_base = IUPAC_dict[call]

        # process "2/2, 1/2, etc."
        if int(allele1) >= 1 and int(allele2) > 1 :

            # deal with homozygotes
            if allele1 == allele2:
                called_base = alt.split(',')[int(allele1) - 1]

            # deal with heterozygotes
            else:

                if IUPAC_ambiguities == False:
                    called_base = 'N'

                else:
                    ref = alt.split(',')[int(allele1) - 1]
                    alt = alt.split(',')[int(allele2) - 1]
                    call = [ref, alt]
                    call.sort()
                    call = tuple(call)
                    called_base = IUPAC_dict[call]

    return called_base

def make_empty_vcf_ordered_dict(vcf_path):
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
    populations_dict = {}
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


def get_slice_indicies(vcf_bgzipped_file, regions, window_size, regions_to_skip=[]):
    """Get slice information from VCF file that is tabix indexed file (bgzipped). """

    # READ IN FILE HEADER
    tbx = pysam.Tabixfile(vcf_bgzipped_file)  # TODO: create try statement to test that file is actually a VCF

    # PARSE LENGTH INFO FROM HEADER

    chrm_lengths = []
    chrm_lengths_dict = {}
    for line in tbx.header:

        if line.startswith("##contig="):

            chrm_name = re.findall(r'ID=.*,', line)
            chrm_name = chrm_name[0].strip('ID=').strip(',')

            chrm_length = re.findall(r'length=.*>', line)
            chrm_length = int(chrm_length[0].strip('length=').strip('>'))

            if chrm_name in regions_to_skip:
                print 'skipping', chrm_name
                continue

            chrm_lengths.append((chrm_name, 1, chrm_length))
            chrm_lengths_dict[chrm_name] = chrm_length

    chrm_lengths = tuple(chrm_lengths)
    tbx.close()

    # GENERATE SLICES
    # current does not make overlapping slices.

    if regions == None:

        for chrm, start, stop in chrm_lengths:

            slice_indicies = itertools.islice(xrange(start, stop + 1), 0, stop + 1, window_size)

            for count, s in enumerate(pairwise(slice_indicies)):
                yield (chrm, s[0], s[1] - 1)  # subtract one to prevent 1 bp overlap
    else:
        for r in regions:
            split_region = re.split(r':|-', r)

            if len(split_region) == 3:
                chrm, start, stop = split_region
            else:
                chrm = split_region[0]
                start, stop = 1, chrm_lengths_dict[chrm]

            start, stop = int(start), int(stop)

            slice_indicies = itertools.islice(xrange(start, stop + 1), 0, stop + 1, window_size)
            for count, s in enumerate(pairwise(slice_indicies)):
                yield (chrm, s[0], s[1] - 1)   # subtract one to prevent 1 bp overlap


def slice_vcf(vcf_bgzipped_file, chrm, start, stop):

    tbx = pysam.Tabixfile(vcf_bgzipped_file)

    try:
        vcf_slice = tbx.fetch(chrm, start, stop)
    except ValueError:
        return None

    else:
        return tuple(row for row in vcf_slice)


def parse_info_field(info_field):

    info_dict = {}
    for item in info_field.split(';'):
        pair = item.split("=")
        if len(pair) == 2:
            info_dict[pair[0]] = pair[1]   # this could be improved on
    return info_dict


def parse_vcf_line(pos, vcf_line_dict):
    """Read in VCF line and convert it to an OrderedDict"""

    pos_parts = pos.strip().split()

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


def get_population_sizes(vcfline, populations):

    sample_counts = {}

    for pop in populations.keys():
        sample_count = 0

        for sample in populations[pop]:
            if vcfline[sample] is not None:
                sample_count += 1

        sample_counts[pop] = sample_count

    return sample_counts


def summarize_population_sizes(dict_of_sizes):
    results = {}
    for pop, sizes, in dict_of_sizes.iteritems():
        sizes = np.array(sizes)
        results[pop + '.sample_count.mean'] = sizes.mean()
        results[pop + '.sample_count.stdev'] = np.std(sizes)

    return results


def pop_size_statistics_2_sorted_list(pop_size_statistics, order):

    if len(order) == 0:
        order = pop_size_statistics.keys()
        order.sort()

    stats = []
    for key in order:
        stat = pop_size_statistics[key]
        stats.append(stat)

    return (stats, order)


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

    #allele_counts = defaultdict({0:0.0,1:0.0,2:0.0,3:0.0,4:0.0})
    allele_counts = dict((key, {0:0.0, 1:0.0, 2:0.0, 3:0.0}) for key in populations.keys())

    for population in populations.keys():
        for sample_id in populations[population]:

            if vcf_line_dict[sample_id] != None:

                genotype = vcf_line_dict[sample_id]
                genotype = genotype["GT"].split("/")   # TODO add phased logic

                if genotype == [".", "."]:
                    continue

                genotype = [int(item) for item in genotype]

                for allele in genotype:
                    allele_counts[population][allele] += 1.0

    return allele_counts

def calc_fstats(allele_counts):

    # CALCULATE ALLELE FREQUENCIES
    #allele_freqs_dict = populations.fromkeys(populations.keys(), None)
    allele_freqs_dict = defaultdict(None)
    populations = allele_counts.keys()

    for population in allele_counts.keys():
        counts = allele_counts[population].values()

        if sum(counts) == 0.0:
            freqs = [0.0] * 4

        else:
            freqs = [count / sum(counts) for count in counts]

        allele_freqs_dict[population] = freqs

    allele_freqs = allele_freqs_dict.values()


    # CACULATE PAIRWISE F-STATISTICS
    pairwise_results = {}

    for population_pair in combinations(populations, 2):

        pop1, pop2 = population_pair
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

            # FORMAT
            values = [Hs_est_, Ht_est_, Gst_est_, G_prime_st_est_,
                      G_double_prime_st_est_, D_est_]
            values_dict = dict(zip(['Hs_est', 'Ht_est', 'Gst_est', \
                                    'G_prime_st_est', 'G_double_prime_st_est',
                                    'D_est'],
                                    values))

            pairwise_results[population_pair] = values_dict

    return pairwise_results


def calc_multilocus_f_statistics(Hs_est_dict, Ht_est_dict):

    multilocus_f_statistics = {}

    for key in Hs_est_dict.keys():

        Hs_est_list = Hs_est_dict[key]
        Ht_est_list = Ht_est_dict[key]

        pairs = zip(Hs_est_list, Ht_est_list)
        pairs = [pair for pair in pairs if np.nan not in pair]

        multilocus_f_statistics[key] = None

        if len(pairs) != 0:

            # THIS REMOVES NaNs FROM THE Ht and Hs LISTS
            Ht_est_list_no_NaN = fstats.de_NaN_list(Ht_est_list)
            Hs_est_list_no_NaN = fstats.de_NaN_list(Hs_est_list)

            n = 2 # fix this, assumes pairs = paired population
            Gst_est = fstats.multilocus_Gst_est(Ht_est_list_no_NaN, Hs_est_list_no_NaN)
            G_prime_st_est = fstats.multilocus_G_prime_st_est(Ht_est_list_no_NaN, Hs_est_list_no_NaN, n)
            G_double_prime_st_est = fstats.multilocus_G_double_prime_st_est(Ht_est_list_no_NaN, Hs_est_list_no_NaN, n)

            # NOTE that the Dest calculation handles NaN's better and
            # can use the raw Hs and Ht calls
            D_est = fstats.multilocus_D_est(Ht_est_list, Hs_est_list, n)

            values_dict = dict(zip(['Gst_est', 'G_prime_st_est', 'G_double_prime_st_est', 'D_est'],\
                          [Gst_est, G_prime_st_est, G_double_prime_st_est, D_est]))

            multilocus_f_statistics[key] = values_dict

    return multilocus_f_statistics


def update_Hs_and_Ht_dicts(f_statistics, Hs_est_dict, Ht_est_dict):

    for pop_pair in f_statistics.keys():

        if Hs_est_dict.has_key(pop_pair) == True:
            Hs_est_dict[pop_pair].append(f_statistics[pop_pair]['Hs_est'])
        else:
            Hs_est_dict[pop_pair] = [f_statistics[pop_pair]['Hs_est']]

        if Ht_est_dict.has_key(pop_pair) == True:
            Ht_est_dict[pop_pair].append(f_statistics[pop_pair]['Ht_est'])
        else:
            Ht_est_dict[pop_pair] = [f_statistics[pop_pair]['Ht_est']]

    return (Hs_est_dict, Ht_est_dict)

def f_statistics_2_sorted_list(multilocus_f_statistics, order=[]):

    if len(order) == 0:

        for pair in multilocus_f_statistics.keys():
            joined_pair = '.'.join(pair)
            for stat in multilocus_f_statistics[pair].keys():
                order.append('.'.join((joined_pair, stat)))

        order.sort()


    stats = []
    for key in order:
        pop1, pop2, stat = key.split(".")

        if multilocus_f_statistics[(pop1, pop2)] == None:                # TO DO: figure out why this occurs
            stats.append(np.nan)
        else:
            value = multilocus_f_statistics[(pop1, pop2)][stat]
            stats.append(value)

    return (stats, order)

def process_header(tabix_file):

    chrm_lenghts_dict = {}
    tabix_file = pysam.Tabixfile(tabix_file)
    
    for line in tabix_file.header:
        
        if line.startswith("##contig") == True: 
            chrm, length = re.split(r"<ID=|,length=", line)[1:]
            length =  int(length.strip(">"))
            chrm_lenghts_dict[chrm] = length

    return chrm_lenghts_dict


def generate_fstats_from_vcf_slices(slice_indicies, populations, header, args):

    for count, si in enumerate(slice_indicies):
        chrm, start, stop = si

        yield [slice_vcf(args.input, chrm, start, stop),
               chrm, start, stop, populations, header,
               args.min_samples]


def process_outgroup(vcf_line, populations):
    og = populations["Outgroup"]

    # Create list of outgroup genotypes
    gt = []
    for s in og:
        if vcf_line[s] != None:
            gt.append(vcf_line[s]["GT"])

    # Get unique genotypes and only return
    # genotypic information if the genotype is
    # homozygous for both samples with

    gt = set(gt)
    if len(gt) == 1:
        gt = tuple(gt)[0]
        gt = set(re.split(r"/|\|", gt))   # split on "/" or "|"

        if len(gt) == 1:
            return tuple(gt)[0]
        else:
            return None

    else:
        return None


def identify_fixed_populations(allele_counts, order):

    if len(order) == 0:
        order = allele_counts.keys()
        order.sort()

    fixed_pops = []
    for pop in order:
        a_values = allele_counts[pop].values()
        non_zero_a = set([a for a in a_values if a != 0])

        if len(non_zero_a) == 1:
            fixed_pops.append(1)
        else:
            fixed_pops.append(0)

    return (fixed_pops, order)


def vcf_line_to_snp_array(vcf_line_dict):
    pass


def calc_slice_stats(data):
    """Main function for caculating statistics.

       Make it easy to add more statistics.
    """

    tabix_slice, chrm, start, stop, populations, header, min_samples = data

    #progress_meter(starting_time, chrm, stop, bp_processed, total_bp_in_dataset)

    if tabix_slice == None or len(tabix_slice) == 0:  # skip empty alignments
        return None

    else:
        # Create lists to store values for final printing.
        output_line = [chrm, start, stop]
        total_depth = []
        snp_wise_f_statistics = []
        population_sizes = defaultdict(list)

        Hs_est_dict = {}
        Ht_est_dict = {}
        snp_count = 0

        for count, line in enumerate(tabix_slice):

            vcf_line_dict = parse_vcf_line(line, header)

            # CREATE FILTERS HERE:
            if vcf_line_dict["FILTER"] != 'PASS':
                continue

            # COUNT SAMPLES IN EACH POPULATION
            for pop, size in get_population_sizes(vcf_line_dict, populations).iteritems():
                population_sizes[pop].append(size)

            # CALCULATE SNPWISE F-STATISTICS
            allele_counts = calc_allele_counts(populations, vcf_line_dict)
            f_statistics = calc_fstats(allele_counts)

            # UPDATE Hs AND Ht DICTIIONARIES
            Hs_est_dict, Ht_est_dict = update_Hs_and_Ht_dicts(f_statistics, Hs_est_dict, Ht_est_dict)
            f_statistics['LOCATION'] = (chrm, start, stop)
            snp_wise_f_statistics.append(f_statistics)

            total_depth.append(int(vcf_line_dict['INFO']["DP"]))
            snp_count = count

        # SUMMARIZE POPULATION WIDE STATISTICS
        pop_size_statistics = summarize_population_sizes(population_sizes)
        multilocus_f_statistics = calc_multilocus_f_statistics(Hs_est_dict, Ht_est_dict)

        # SKIP SAMPLES WITH TOO MANY NANs

        if len(multilocus_f_statistics.values()) == 0:
            return None

        elif multilocus_f_statistics.values()[0] is None:
            return None

        else:
            # UPDATE OUTPUT LINE WITH DEPTH INFO
            total_depth = np.array(total_depth)
            output_line += [snp_count, total_depth.mean(), total_depth.std()]

            return ([chrm, start, stop, snp_count, total_depth.mean(), total_depth.std()], \
                 pop_size_statistics, multilocus_f_statistics)

