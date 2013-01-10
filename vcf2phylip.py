#!/usr/bin/env python
# encoding: utf-8
"""
vcf2phylip.py

Created by Nick Crawford on 2011-11-18.
Copyright (c) 2011

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses

The author may be contacted at ngcrawford@gmail.com
"""

import os
import re
import sys
import VCF
import gzip
import pysam
import shlex
import random
import argparse
import itertools
import numpy as np
from subprocess import Popen, PIPE
from collections import namedtuple

def get_args():
    """Parse sys.argv"""
    parser = argparse.ArgumentParser()
    parser.add_argument('-r','--regions',
                        required=True,
                        type=str,
                        help="Chromosomal region in the format: 'Chrm:start-stop'")
                        
    parser.add_argument('-o','--output',
                        type=argparse.FileType('w'),
                        default=sys.stdout,
                        help='Path to output. (default is STOUT)')

    parser.add_argument('-b','--bootstraps', 
                        type=int,
                        help='Calculate bootstraps.')

    parser.add_argument('input', 
                        nargs=1, 
                        help='bgzipped and indexed VCF file')
    
    args = parser.parse_args()
    return args

def makeDataTuple(vcf):
    """Setup a labeled tuple to store the data."""

    chrm_data = {}

    for count, line in enumerate(vcf):

        if "##contig" in line:
            contigline = line.split("<")[1]
            contigline = contigline.strip(">\n").split(",")
            contigline = [item.split("=")[1] for item in contigline]
            chrm_data[contigline[0]] = int(contigline[1])

        if line.startswith("#CHROM"):
            field_labels = line.strip().strip("#").split("\t")
            field_labels = [item.strip().split('.')[0].replace("-", "_") for item in field_labels]
            break

    position_data = namedtuple('base', field_labels)
    return (position_data, chrm_data)

def array2OnelinerAlignment(info, taxa, bases):
    """Convert array of array of taxa and an array of bases to one-liner."""

    oneliner = info
    for count, seq in enumerate(bases):
        oneliner += taxa[count]+","+''.join(itertools.chain(bases[count])) + ","
    oneliner = oneliner[:-1] + ";"
    return oneliner

def callSNPs(current_base, numb_of_seqs):
    """Call the SNPs. Duh!"""

    blanks =  np.zeros(numb_of_seqs, np.string0)

    if current_base.FILTER == 'LowQual':
        blanks.fill("-")

    if current_base.FORMAT == 'GT':
        blanks.fill("-")

    for count, snp_call in enumerate(current_base[9:]):
        base = VCF.process_snp_call(snp_call, current_base.REF, current_base.ALT)
        blanks[count] = base

    return blanks

def count_informative_sites(alignment_array):
    """Informative Sites must have two different SNPs"""
    informative_sites = 0
    for site in alignment_array:
        unique_sites = set(site)
        if len(unique_sites) >= 3:
            informative_sites += 1
    return informative_sites

def get_subset_vcf(chrm, start, stop):
    base_dir = "/Users/MullenLab/Desktop/Grad_Students/Nick/butterfly_practice"
    cli = """java -Xmx6g -jar /Users/MullenLab/Source/gatk/dist/GenomeAnalysisTK.jar \
      -R {3}/Butterfly.merge.scaffolds.fa  \
      -T SelectVariants \
      --variant {3}/32_butterflies_vcfs/32.butterflies.good_BGI_snps.combined.vcf \
      -L {0}:{1}-{2}""".format(chrm, start, stop, base_dir)

    cli_parts = shlex.split(cli)
    vcf = Popen(cli_parts, stdin=PIPE, stderr=PIPE, stdout=PIPE).communicate()[0]
    return vcf

def generate_bootstraps(chrm, chrm_len, window_size, numb_of_reps):

    start_sites = [random.randint(0, chrm_len-window_size) for item in range(numb_of_reps)]
    replicates = []
    for count, start in enumerate(start_sites):
        vcf = get_subset_vcf(chrm, start, start+window_size)
        yield vcf

def parse_window_vcf(vcf, start, stop, window_size, chrm, fout):
    # SETUP NAMED TUPLE TO STORE INFO FROM A SINGLE BASE
    field_labels = []
    position_data, chrm_data = makeDataTuple(vcf)

    # SETUP MULTIPLE ALIGNMENT ARRAY
    numb_of_seqs = len(position_data._fields[9:])
    alignment = np.zeros((window_size,numb_of_seqs), np.string0)

    # SETUP COUNTERS
    current_base = None
    current_window = 1
    line_count = 0
    windows = range(0, chrm_data[chrm], window_size)
    current_data = []
    informative_sites = []

    # PARSE VCF FIlE
    snp_count = 0

    for line in vcf:
        # SKIP HEADER
        if line.startswith("#CHROM"):
            line_count = 0

        # START PROCESSING ALIGNED BASES
        if line.startswith(chrm):
            current_base = position_data._make(line.strip().split("\t"))
            base_calls = callSNPs(current_base, numb_of_seqs)
            current_data.append(base_calls.copy())

    alignment = np.array(current_data)
    inform_sites = count_informative_sites(alignment)
    if current_base == None:
        return 'error'
    else:
        taxa = current_base._fields[9:]
        info = 'chrm={0},start={1},stop={2},inform_sites={3}'.format(current_base.CHROM, start, stop, inform_sites)
        oneliner = array2OnelinerAlignment(info, taxa, alignment.transpose())

    if ":" in oneliner and oneliner[-1] == ';': # this prevents bad alignments from getting printed
        return oneliner
    else:
        return 'error'

def header_slices(vcf, window_size=5000):
    vcf = pysam.Tabixfile(vcf)
    
    slices = {}
    for line in vcf.header:
        
        if line.startswith('##contig'):
            
            line_dict = dict([item.split("=") for item in line.split("<")[1].strip('>').split(",")])
            length = int(line_dict['length'])
            
            if length < window_size: continue

            start = (length % window_size)/2
            stop = ((length/window_size) * window_size) + start
            
            s = xrange(start,stop,window_size)
            slices[line_dict['ID']] = s
    
    return slices


def slice_vcf(vcf, chrm, start, stop):
    cli = """tabix {0} {1}:{2}-{3}""".format(vcf, chrm, start, stop)
    cli_parts = shlex.split(cli)
    vcf = Popen(cli_parts, stdin=PIPE, stderr=PIPE, stdout=PIPE).communicate()[0]
    return vcf

def process_vcf_slice(tabix_file, chrm, start, stop, position_data):

    tbx = pysam.Tabixfile(tabix_file)
    tbx_lines = tbx.fetch(chrm, start, stop)

    numb_of_seqs = len(position_data._fields[9:])
    alignment = np.zeros((stop-start,numb_of_seqs), np.string0)

    # This 'error handling' needs to be rewritten.
    current_data = []
    if tbx_lines == None:
        return 'error'

    for line in tbx_lines:
        current_base = position_data._make(line.strip().split("\t"))
        base_calls = callSNPs(current_base, numb_of_seqs)
        current_data.append(base_calls.copy())

    alignment = np.array(current_data)
    inform_sites = count_informative_sites(alignment)
    
    if current_base == None:
        return 'error'
    else:
        taxa = current_base._fields[9:]
        info = "tree 'chrm={0},start={1},stop={2},inform_sites={3}':".format(current_base.CHROM, start, stop, inform_sites)
        oneliner = array2OnelinerAlignment(info, taxa, alignment.transpose())

    if ":" in oneliner and oneliner[-1] == ';': # this prevents bad alignments from getting printed
        return oneliner
    else:
        return 'error'

def oneliner2phylip(line):
    """Convert one-liner to phylip format."""
    seqs = line.strip(";").split(':')[-1].split(',')
    label_seqs = zip(seqs[:-1:2],seqs[1::2])
    taxa_count = len(label_seqs)
    seq_length = len(label_seqs[0][1])
    alignment = "%s %s\n" % (taxa_count, seq_length) # add header
    for taxa_name, seq in label_seqs:
        taxa_name = taxa_name.strip()
        alignment += '%-10s%s\n' % (taxa_name, seq)
    return alignment

def main():

    # SETUP ARGS
    args = get_args()

    chrm, start, stop = re.split(r':|-', args.regions)
    start, stop = int(start), int(stop)
    
    # OPEN VCF
    vcf_file = gzip.open(args.input[0],'rb')
    position_data, chrm_data = makeDataTuple(vcf_file)
    vcf_file.close()

    oneliner = process_vcf_slice(args.input[0], chrm, start, stop, position_data)
    oneliner = oneliner2phylip(oneliner)
    if oneliner != 'error':
        sys.stdout.write(oneliner+"\n")

if __name__ == '__main__':
    main()




