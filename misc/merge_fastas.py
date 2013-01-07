#!/usr/bin/env python
# encoding: utf-8

import os
import sys
import argparse
try:
    from Bio import SeqIO
except:
    sys.write('Bio Python needs to be installed: biopython.org\n')
    sys.exit(0)

def default_args():
    """Parse sys.argv"""
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input', 
                        required=True, 
                        type=str,
                        help='Path to fasta file.')

    parser.add_argument('-n', '--number-of-Ns', 
                        default=300, 
                        type=int,
                        help="""Number of Ns to append between sequences.""")

    parser = parser.parse_args()
    return parser

def process_output():
    pass


args = default_args()

# Create file with output index positions
fout_idx = os.path.splitext(args.input)[0] + ".merged_contigs.positions.txt"
fout = os.path.splitext(args.input)[0] + ".merged_contigs.fa"
Ns = 'N' * args.number_of_Ns

handle = open(args.input, "rU")
for record in SeqIO.parse(handle, "fasta") :

    new_seq = len(record.seq + Ns)
    print new_seq % 60


handle.close()