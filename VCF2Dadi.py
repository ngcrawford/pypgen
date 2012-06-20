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


import VCF
import argparse
from copy import copy

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
                        help='Chrm:start-stop')


    args = parser.parse_args()

    populations_dict  = {}
    for pop in args.populations:
        pop_name, sample_ids = pop.strip().split(":")
        sample_ids = sample_ids.split(",")
        populations_dict[pop_name] = sample_ids

    args.populations = populations_dict

    chrm = [args.region.split(":")[0]]
    start_stop = [int(item) for item in args.region.split(":")[1].split("-")]
    args.region = chrm + start_stop

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
	allele_counts = row["outgroups"]
	
	if allele_counts['REF'] != 0:
		outgroup_base = raw_calls['REF']
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


def main(args):

	vcf = VCF.VCF()
	vcf.populations = args.populations
	vcf.set_header(args.input)
	chunk = vcf.vcf_chunk_2_dadi(args.input, args.populations, *args.region)

	g = count_alleles(chunk, args.populations)

	pop_ids = args.populations.keys()[:-1]

	# TO DO write function to create header
	dadi_header = ' '.join(['Outgroup','Helis','Allele1',pop_ids[0],pop_ids[1],pop_ids[2],\
	       'Allele2',pop_ids[0],pop_ids[1],pop_ids[2],'Chrm','Pos'])


	fout = open(args.output,'w')
	fout.write(dadi_header + "\n")

	for row_count, row in enumerate(g):
		raw_calls = chunk[row_count] 

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

if __name__ == '__main__':
	args = get_args()
	main(args)


