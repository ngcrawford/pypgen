#!/usr/bin/env python
# encoding: utf-8

import re
import textwrap
import multiprocessing
from pypgen.parser.VCF import *
from pypgen.misc.helpers import *


def process_header(tabix_file):

    chrm_lenghts_dict = {}
    tabix_file = pysam.Tabixfile(tabix_file)

    for line in tabix_file.header:

        if line.startswith("##contig") == True:
            chrm, length = re.split(r"<ID=|,length=", line)[1:]
            length = int(length.strip(">"))
            chrm_lenghts_dict[chrm] = length

    return chrm_lenghts_dict


def main():
    # get args.
    args = default_args()


    args.description  = textwrap.dedent("""\
    vcf_snpwise_fstats.py version 0.2.0 beta by 
    Nicholas Crawford (ngcrawford@gmail.com)

    Working Example:
        python vcf_snpwise_fstats.py \\
        -i data/example.vcf.gz \\
        -p outgroups:h665,i02-210 \\
        pop1:c511,c512,c513,c514,c515,c563,c614,c630,c639,c640 \\
        pop2:m523,m524,m525,m589,m675,m676,m682,m683,m687,m689 \\
        -c 2 \\
        -r Chr01:1-10001 | head
    """)

    args = args.parse_args()


    # TODO:
    # test that pysam is installed.
    # bgzip check. MDSum?

    # Calculate the total size of the dataset
    chrm_lengths = process_header(args.input)

    # Get information about samples from the header.
    # this becomes the precursor to the VCF row
    empty_vcf_line = make_empty_vcf_ordered_dict(args.input)

    # Convert populations input into a dict of pops where
    # values are lists of samples
    populations = parse_populations_list(args.populations)

    fstat_order = []    # store order of paired samples.
    pop_size_order = []
    fixed_alleles_order = []
    vcf_count = 0

    for count, line in enumerate(open_vcf(args)):

        if line.startswith('#') == True:
            continue

        vcf_line = parse_vcf_line(line, empty_vcf_line)

        # APPLY BASIC FILTER
        if vcf_line["FILTER"] != "PASS":
            continue

        # SKIP ALLELES WITH NO INFORMATION
        # NOTE: with multi alleleic snps this may
        # still leave some populations with nans..
        # but should account for most issues
        if ',' not in vcf_line["INFO"]["AF"]:
            if float(vcf_line["INFO"]["AF"]) == 1.0:
                continue

        vcf_count += 1

        allele_counts = calc_allele_counts(populations, vcf_line)

        fixed_alleles, fixed_alleles_order = identify_fixed_populations(allele_counts, fixed_alleles_order)
        fstats = calc_fstats(allele_counts)
        pop_size_stats = get_population_sizes(vcf_line, populations)

        f_stats, fstat_order = f_statistics_2_sorted_list(fstats, order=fstat_order)
        pop_size_stats, pop_size_order = pop_size_statistics_2_sorted_list(pop_size_stats, order=pop_size_order)

        chrm = vcf_line['CHROM']
        pos = vcf_line['POS']

        pop_size_stats = [float_2_string(i, places=4) for i in pop_size_stats]
        f_stats = [float_2_string(i, places=4) for i in f_stats]

        if vcf_count == 1:
            args.output.write(','.join(['chrm', 'pos'] \
                       + map(str, pop_size_order + fstat_order + [pop + "_fixed" for pop in fixed_alleles_order])) + "\n")
            args.output.write(','.join(map(str, [chrm, pos] + pop_size_stats + f_stats + fixed_alleles)) + "\n")

        else:
            args.output.write(','.join(map(str, [chrm, pos] + pop_size_stats + f_stats + fixed_alleles)) + "\n")


        #previous_update_time = progress_meter(previous_update_time, chrm, stop, bp_processed, total_bp_in_dataset)


if __name__ == '__main__':
    main()
