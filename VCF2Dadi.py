#!/usr/bin/env python
# encoding: utf-8

"""
run_dadi.py

Created by Nick Crawford on 2012-07-05.
Copyright (c) 2012

The author may be contacted at ngcrawford@gmail.com


python VCF2Dadi.py \
-i test_data/butterfly.vcf.gz \
-o test_data/butterfly.dadi.input.txt \
-l Chr01:500-1000 \
-p cydno:c511,c512,c513,c514,c515,c563,c614,c630,c639,c640 \
outgroups:h665,i02-210 \
melpo:m523,m524,m525,m589,m675,m676,m682,m683,m687,m689 \
pachi:p516,p517,p518,p519,p520,p591,p596,p690,p694,p696 
"""

import sys
import VCF
import dadi
import argparse
from copy import copy, deepcopy

def get_args():
	"""Parse sys.argv"""
	parser = argparse.ArgumentParser()

	parser.add_argument('-i','--input', required=True,
						help='Path to VCF file.')
	
	parser.add_argument('-o','--output',
						help='Path to output csv file. If path is not set defaults to STDOUT.')

	parser.add_argument('-p','--populations', nargs='+',
						help='Names of populations and samples. The format is: "PopName:sample1,sample2,sample3,etc..."')

	parser.add_argument('-L','--region',default=None, type=str,
						help='chrm:start-stop')

	parser.add_argument('-w', '--window-size', type=int,
						help='The size of the windows')

	parser.add_argument('-overlap', type=int,
						help='The number of base pairs each window overlaps the previous')

	args = parser.parse_args()

	populations_dict  = {}
	for pop in args.populations:
		pop_name, sample_ids = pop.strip().split(":")
		sample_ids = sample_ids.split(",")
		populations_dict[pop_name] = sample_ids

	args.populations = populations_dict

	if args.region != None:
		if len(args.region.split(":")) == 2:
			chrm = [args.region.split(":")[0]]
			start_stop = [int(item) for item in args.region.split(":")[1].split("-")]
			args.region = chrm + start_stop

	else:
		args.region = [args.region]


	return args

def process_snp(snp_call):
	if snp_call == "0/0":
		return (2,0)
	
	elif snp_call == "1/1":
		return (0,2)
	
	elif snp_call == '1/0' or \
		   snp_call == '0/1':
		return (1,1)
	
	# skip multiallelic sites
	else:
		return (0,0)

def count_alleles(chunk, pops):

	results = []
	for line in chunk:
		pop_counts = {}
		for pop in pops.keys():
			allele_counts = {'REF':0, 'ALT':0}
			for sample in pops[pop]:
				if line[sample] != None:
					ref, alt = process_snp(line[sample]['GT'])
					allele_counts['REF'] += ref
					allele_counts['ALT'] += alt

			pop_counts[pop] =  allele_counts.copy()

		results.append(copy(pop_counts))
	
	return results

def make_triplet(base):
	return "-{0}-".format(base)

def check_outgroup(row):
	outgroup = row["outgroups"]
	counts = set([value for value in outgroup.values() if value != 0])
	if len(counts) > 1:
		return False
	else:
		return True

def get_outgroup_base(row, raw_calls):
	
	if row["outgroups"]['REF'] > row["outgroups"]['ALT']:
		outgroup_base = raw_calls['REF']
	
	elif sum(row['outgroups'].values()) == 0 : # e.g, == {'ALT': 0, 'REF': 0}
		ref = sum([item['REF'] for item in row.values()])
		alt = sum([item['ALT'] for item in row.values()])
		if ref > alt:
			outgroup_base = raw_calls['REF']
		else: 
			outgroup_base = raw_calls['ALT']
	
	else:
		outgroup_base = raw_calls['ALT']	
	
	return outgroup_base

def get_ingroup_major_allele(row, raw_calls, outgroup_allele):
	
	pops = row.keys()
	pops.remove('outgroups')

	ref_sum = sum([row[pop]['REF'] for pop in pops])
	alt_sum = sum([row[pop]['ALT'] for pop in pops])

	if ref_sum < alt_sum:
		return raw_calls['ALT']

	else:
		return raw_calls['REF']		

def calc_in_group_ref(row):
	for key in row.keys().remove('outgroups'):
		print key

def create_dadi_header(args):
	pop_ids = args.populations.keys()
	dadi_header = ['Outgroup','Ingroup','Allele1','Allele2','Chrm','Pos']
	dadi_header[3:3] = pop_ids
	dadi_header[-2:2] = pop_ids
	dadi_header = ' '.join(dadi_header)
	return dadi_header

def make_dadi_fs(args, region):

	vcf = VCF.VCF()
	vcf.populations = args.populations
	vcf.set_header(args.input)

	pop_ids = args.populations.keys()

	# Get slice and setup output dictionaries
	chunk = vcf.slice_vcf(args.input, *region)
	if chunk == None:
		return None

	else:
		g = count_alleles(chunk, args.populations)

		final_dadi = {}
		population_level_dadis = dict.fromkeys(pop_ids,{})

		for row_count, row in enumerate(g):

			raw_calls = chunk[row_count] 
			row['outgroups'] = {'ALT': 0, 'REF': 0} # set empty outgroup

			# To Do: Need to create a function to fill outgroup if one is defined.
			# The heliconius dataset, for example, has this.

			if check_outgroup(row) == False: continue # skip outgroup not fixed at one value
			if len(raw_calls['REF']) > 1 or len(raw_calls["ALT"]) > 1: continue # skip multi allelic sites
			
			# CALL BASE FOR OUTGROUP
			outgroup_allele = get_outgroup_base(row, raw_calls)

			# CALL MAJOR ALLELE (BASE) FOR INGROUP
			major_allele = get_ingroup_major_allele(row, raw_calls, outgroup_allele)

			# POLORIZE REF AND ALT FOR INGROUP
			if major_allele != raw_calls['REF']:
				ref, alt = ('ALT','REF')
			else:
				ref, alt = ('REF','ALT')


			calls = {}
			for count, pop in enumerate(pop_ids):
				calls[pop] = (row[pop][ref], row[pop][alt])

			row_id = "{0}_{1}".format(raw_calls['CHROM'],raw_calls['POS'])
		
			dadi_site = {'calls': calls,
				   'context': make_triplet(major_allele),
				   'outgroup_context': make_triplet(outgroup_allele),
				   'outgroup_allele': outgroup_allele,
				   'segregating': (raw_calls[ref], raw_calls[alt])
				   }

			final_dadi[row_id] = dadi_site

		return (final_dadi, pop_ids)

def generate_slices(args):

	vcf = VCF.VCF()
	vcf.populations = args.populations
	vcf.set_chrms(args.input)

	chrm_2_windows = vcf.chrm2length.fromkeys(vcf.chrm2length.keys(),None)
	
	for count, chrm in enumerate(vcf.chrm2length.keys()):

		length = vcf.chrm2length[chrm]
		window_size = args.window_size
		overlap = args.overlap

		# Skip contigs that are to short
		if length <= window_size: continue
		
		# Fit windows into remaining space
		if (length % window_size) > overlap:
			start = (length % window_size)/2
			stop = (length - window_size) - overlap/2

		# Prevent windows from invading remaining space 
		if (length % window_size) <= overlap:
			start = (length % window_size)/2
			stop = length - overlap*2
				
		starts = range(start, stop, overlap)
		stops = [i+window_size for i in starts]
		windows = zip(starts, stops)
		
		chrm_2_windows[chrm] = windows

	return chrm_2_windows

def main(args):

	vcf = VCF.VCF()
	vcf.populations = args.populations
	vcf.set_header(args.input)

	pop_ids = args.populations.keys()

	# get slice and setup output dictionaries
	chunk = vcf.vcf_chunk_2_dadi(args.input, args.populations, *args.region)
	g = count_alleles(chunk, args.populations)


	# Create Header Row
	dadi_header = create_dadi_header(args)

	fout = open(args.output,'w')
	fout.write(dadi_header + "\n")

	for row_count, row in enumerate(g):

		raw_calls = chunk[row_count] 
		row['outgroups'] = {'ALT': 0, 'REF': 0} # set empty outgroup

		# To Do: Need to create a function to fill outgroup if one is defined.
		# The heliconius dataset, for example, has this.

		if check_outgroup(row) == False: continue # skip outgroup not fixed at one value
		if len(raw_calls['REF']) > 1 or len(raw_calls["ALT"]) > 1: continue # skip multi allelic sites
		
		# CALL BASE FOR OUTGROUP
		outgroup_allele = get_outgroup_base(row, raw_calls)

		# CALL MAJOR ALLELE (BASE) FOR INGROUP
		major_allele = get_ingroup_major_allele(row, raw_calls, outgroup_allele)

		# POLORIZE REF AND ALT FOR INGROUP
		if major_allele != raw_calls['REF']:
			ref, alt = ('ALT','REF')
		else:
			ref, alt = ('REF','ALT')

		#  CREATE DADI ROW
		dadi_row = [make_triplet(major_allele), make_triplet(outgroup_allele)]

		for count, pop in enumerate(pop_ids):
			if count == 0:
				dadi_row.append(chunk[row_count][ref])
				dadi_row.append(row[pop][ref])
			else:
				dadi_row.append(row[pop][ref])


		for count, pop in enumerate(pop_ids):
			if count == 0:
				dadi_row.append(chunk[row_count][alt]) 
				dadi_row.append(row[pop][alt])
			else:
				dadi_row.append(row[pop][alt])

		dadi_row.append(raw_calls['CHROM'])
		dadi_row.append(raw_calls['POS'])

		dadi_row = " ".join([str(item) for item in dadi_row])
		fout.write(dadi_row + "\n")

def create_header(pop_ids):

	pop_values = ['Tajimas_D','W_theta','pi','Seg_Sites']

	final_header = ['chrm', 'start','stop','Fst']
	for count, pop in enumerate(pop_ids):
		final_header += [pop+ "." + i for i in pop_values]

	return final_header


def sliding_window_dadi(args):

	# Using the data in the VCF header generate all slices
	print 'Generating Slices...'
	slices = generate_slices(args)
	
	# Open output file
	print 'Processing Slices...'
	fout = open(args.output,'w')

	for key_count, chrm in enumerate(slices.keys()):

		if args.region != [None]:
			chrm = args.region[0]
			#if key_count == 2: break
		
		if slices[chrm] == None:
			print "Skipped:", chrm, ' - no slices'
			continue

		for count, s in enumerate(slices[chrm]):

			# Break out of loop if loop proceeds beyond
			#   defined region (-L=1:1-5000 = 5000)
			if s[-1] > args.region[-1] and args.region[-1] != None: break

			region = [chrm] + list(s)

			# Setup Pairwise Dadi
			dadi_data = make_dadi_fs(args, region)
			if dadi_data == None: continue # skip empty calls

			dd, pop_ids = dadi_data
			projection_size = 10
			pairwise_fs  = dadi.Spectrum.from_data_dict(dd, pop_ids, [projection_size]*2)

			# Write output header
			if count == 0: fout.write(','.join(create_header(pop_ids)) + "\n")

			# Create final line, add Fst info
			final_line = region
			final_line += [pairwise_fs.Fst()]

			# Add in population level stats
			for pop in pop_ids:
				fs = dadi.Spectrum.from_data_dict(dd, [pop], [projection_size])
				final_line += [fs.Tajima_D(), fs.Watterson_theta(), fs.pi(), fs.S()]

			# write output
			final_line = [str(i) for i in final_line]
			fout.write(','.join(final_line) + "\n")

		# Don't process any more keys than necessary
		print 'Processed {0} of {1} slices from contig {2}'.format(count+1, len(slices[chrm]), chrm)
		if args.region != [None]: break

	fout.close()


# def test_slices(args):
# 	vcf = VCF.VCF()
# 	vcf.set_header(args.input)
# 	print 'testing slices'
# 	slices = generate_slices(args)
# 	input_file = "/usr3/graduate/ngcrawfo/genomics/bulk-seg-anoMar/vcf/CAP_MAR.di_allelic.stampy.vcf.gz"
# 	for chrm in slices.keys():
# 		print chrm, len(vcf.slice_vcf(input_file, chrm, start=1, stop=2500))


if __name__ == '__main__':
	args = get_args()
	#test_slices(args)
	sliding_window_dadi(args)





