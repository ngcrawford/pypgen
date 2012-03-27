#!/usr/local/epd7/bin/python python
# encoding: utf-8
"""
sync2Fst.py

Created by Nicholas Crawford (c) 2010 Nicholas Crawford. All rights reserved.
"""

import os
import sys
import gzip
import pandas
import argparse
import numpy as np
from fstats import *
import multiprocessing
from textwrap import dedent
from itertools import izip_longest

def interface():
	p = argparse.ArgumentParser()
	p.add_argument('-i','--input-file',
        help='Path to input file.')
	p.add_argument('-o','--output-file',
        help='Path to output file.')
	p.add_argument('-w','--window-size',
		type=int,
		help='Size of sliding window (bp)')
	args = p.parse_args()
 	return args

def no_window():
	
	args = interface()

	Ns = [10,8]
	n = float(len(Ns))
	Ns_harm = harmonic_mean(Ns)
	Ns_harm_chao = harmonic_mean_chao(Ns)

	fin = open(args.input_file,'rU')
	fout =  open(args.output_file,'w')
	
	for count, line in enumerate(fin):
		line_parts = line.strip().split()
		
		# CALCULATE BASIC STATS
		chrm, pos, refbase = line_parts[:3]
		pops = line_parts[3:]
		allele_freqs, coverage = calc_SNP_freqs(pops)
		min_coverage = min(coverage)
		max_coverage = max(coverage)

		# CALCULATE Hs AND Ht
		Hs_prime_est_ = Hs_prime_est(allele_freqs,n)
		Ht_prime_est_ = Ht_prime_est(allele_freqs,n)
		Hs_est_ = Hs_est(Hs_prime_est_, Ns_harm)
		Ht_est_ = Ht_est(Ht_prime_est_,Hs_est_,Ns_harm,n)

		# CALCULATE F-STATISTICS
		Gst_est_ = Gst_est(Ht_est_, Hs_est_)
		G_prime_st_est_ = G_prime_st_est(Ht_est_, Hs_est_, Gst_est_, n)
		G_double_prime_st_est_ = G_double_prime_st_est(Ht_est_, Hs_est_, n)
		D_est_ = D_est(Ht_est_, Hs_est_, n)
		
		# PRINT OUTPUT
		values = [chrm, pos, refbase, min_coverage, max_coverage, Hs_est_, \
		Ht_est_, Gst_est_, G_prime_st_est_, G_double_prime_st_est_, D_est_,]
		result =  ','.join([str(item) for item in values])+"\n"
		fout.write(result)

	fin.close()
	fout.close()



def calc_values(allele_freqs):
	n = 2
	Ns = [allele_freqs[0].sum(), allele_freqs[1].sum()]
	Ns_harm = harmonic_mean(Ns)
	Hs_prime_est_ = Hs_prime_est(allele_freqs,n)
	Ht_prime_est_ = Ht_prime_est(allele_freqs,n)
	Hs_est_ = Hs_est(Hs_prime_est_, Ns_harm)
	Ht_est_ = Ht_est(Ht_prime_est_,Hs_est_,Ns_harm,n)

	# CALCULATE F-STATISTICS
	Gst_est_ = Gst_est(Ht_est_, Hs_est_)
	G_prime_st_est_ = G_prime_st_est(Ht_est_, Hs_est_, Gst_est_, n)
	G_double_prime_st_est_ = G_double_prime_st_est(Ht_est_, Hs_est_, n)
	D_est_ = D_est(Ht_est_, Hs_est_, n)
	
	# PRINT OUTPUT
	return [Hs_prime_est_, Ht_prime_est_, Hs_est_, Ht_est_, Gst_est_, G_prime_st_est_, D_est_ ]


def process_window(window):

	alleles = np.array(window.alleles)
	alleles = np.array([gt for gt in alleles])

	# Calculate Fstats 
	f_stats = np.array([calc_values(pop_pair) for pop_pair in alleles \
		if (pop_pair[0] != pop_pair[1]).all() == False])

	means = [str(mean) for mean in f_stats.mean(axis=0)] 
	stdevs = [str(std) for std in f_stats.std(axis=0)]
	return ','.join(means) + ',' + ','.join(stdevs)

def parse_windows(args):

	pool = multiprocessing.Pool(4)

	# handle gzipped files
	if args.input_file.endswith('gz'):
		fin = gzip.open(args.input_file, 'rb')
	else:
		fin = open(args.input_file,'rU')
	
	fout = open(args.output_file,'w')

	current_pos = None
	current_chm = None
	current_chunk = []

	chunks = []
	chunk_count = 0

	for count, line in enumerate(fin):

		line_parts = line.strip().split()
		
		# CALCULATE BASIC STATS
		chrm, pos, refbase = line_parts[:3]
		pos = int(pos)
		pops = line_parts[3:]
		allele_freqs, coverage = calc_SNP_freqs(pops)
		allele_freqs = np.array(allele_freqs)

		min_coverage = min(coverage)
		max_coverage = max(coverage)

		# INITIALIZE CURRENT CHRM AND POS 
		if current_chm == None:
			current_chm = chrm
			current_pos = pos
		
		current_chunk.append([chrm, pos, refbase, min_coverage,\
		                      max_coverage, np.array(allele_freqs)])

		# TEST IF WINDOW HAS BEEN REACHED
		if pos - current_pos >= args.window_size:
			if pos - current_pos != args.window_size: # handle gaps in alignment 
				current_chm = chrm
				current_pos = pos
				continue
			
			# AT THIS POINT YOU'VE GOT THE WINDOW 
			#print current_chm, current_pos, chrm, pos, pos - current_pos
			
			current_chunk = pandas.DataFrame(current_chunk,\
				columns=('chrm','pos','base','min_coverage', \
				         'max_coverage','alleles'))
			
			chunks.append(current_chunk)
			chunk_count += 1
			current_chm = chrm
			current_pos = pos
			current_chunk = []

		# REINITIALIZE IF CHRM CHANGES
		# Note: this will skip/truncate windows at ends of chrms
		if current_chm != chrm:
			current_chm = chrm
			current_pos = pos
		
		if chunk_count % 1000 == 0 and chunk_count != 0: 
			print count, chunk_count, current_chm, current_pos
			fstats = pool.map(process_window, chunks)
			for count, chunk in enumerate(chunks):
				
				line =  [str(item) for item in chunk.ix[0,:3].values] \
				+ [str(chunk.mean()['min_coverage'])] \
				+ [str(chunk.mean()['max_coverage'])] \
				+ [fstats[count]] 
				
				line =','.join(line)
				fout.write(line+"\n")
			
			chunk_count = 0
			chunks = []
	
	# CLEAN UP
	fin.close()
	fout.close()


if __name__ == '__main__':
	args = interface()
	
	if args.window_size == None:
		no_window()
	
	else:
		parse_windows(args)


