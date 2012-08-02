#!/usr/bin/env python
# encoding: utf-8

"""
removeFaFromFasta.py

Created by Nick Crawford on 2012-31-07.

The author may be contacted at ngcrawford@gmail.com

Requires samtools and pysam to be installed.
http://samtools.sourceforge.net/
http://code.google.com/p/pysam/
"""

import os
import pysam
import argparse

def get_args():
	"""Parse sys.argv"""
	parser = argparse.ArgumentParser()

	parser.add_argument('-fa','--fasta', required=True,
						help='Path to fasta file to filter.')

	parser.add_argument('-o','--output', required=True,
						help='Path to filtered fasta file.')

	parser.add_argument('-l', '--filter-list', required=True, nargs='+',
						help='List of of chrm/contig names to remove.')

	args = parser.parse_args()
	return args


def splitIterator(text, size):
	"Iterator that splits string into list of substrings."
	assert size > 0, "size should be > 0"
	for start in xrange(0, len(text), size):
		yield text[start:start + size]

def main(args):
	# Setup files
	fasta = args.fasta
	filtered_fa = args.output
	faidx = args.fasta + '.fai'

	# Make sure .fai exists
	try:
		os.path.exists(faidx)
	except:
		print "You need to index the fasta file with samtools.\n {:} does not exist.".format(faidx)		
	
	# Filter names for contings
	bad_names = args.filter_list
	chrm_names = (line.strip().split()[0] for line in open(faidx,'rU'))
	filtered_chrm_names = (cn for cn in chrm_names if cn not in bad_names)

	# Write contigs/chrms to output fasta
	fasta = pysam.Fastafile(fasta)
	filtered_fa = open(filtered_fa,'w')
	
	for name in filtered_chrm_names:
	
		print "Processed:", name
		chrm = fasta.fetch(name)
		filtered_fa.write(">" + name + "\n")
	
		# split lines on 80 characters
		[filtered_fa.write(chars + "\n" ) for chars in splitIterator(chrm, 80)]

	# Clean up open filesls
	filtered_fa.close()	

if __name__ == '__main__':
	args = get_args()
	main(args)	





