import VCF
import gzip
import argparse
import pysam
from subprocess import Popen, PIPE
from collections import OrderedDict

def split(str, num):
    return [ str[start:start+num] for start in range(0, len(str), num) ]


fasta = "/Users/ngcrawford/Desktop/Genomes/Anolis_carolinensis.AnoCar2.0.67.dna_rm.toplevel.fa"
fa = pysam.Fastafile(fasta)

vcf = VCF.VCF()

input_vcf = "CAP_MAR.50thousand.vcf.gz"
input_vcf = gzip.open(input_vcf,'rb')
modified_fasta = open('CAP_MAR.50thousand.modified.fa','w')

# SETUP COUNTERS
previous_base = {'CHROM':None,'POS':None,'SEQ':None}
line_count = None
sequence = ''

for line in input_vcf:

    # SKIP HEADER
    if line.startswith("#CHROM"):
        vcf.header = line.strip("#").strip().split()
        vcf.__header_dict__ = OrderedDict([(item,None) for item in vcf.header])
        modified_fasta.write(">" + vcf_line['CHROM'] + "\n")
        line_count = 0

    # START PROCESSING ALIGNED BASES
    if line_count > 0:
        vcf_line = vcf.parse_vcf_line(line)
        
        # Reset when VCF when cursor moves to new chrm 
        # and update fasta with new fasta header
        if previous_base['CHROM'] != vcf_line['CHROM'] and previous_base['CHROM'] != None: 
        	previous_base = {'CHROM':None,'POS':None,'SEQ':None}
        	modified_fasta.write(">" + vcf_line['CHROM'] + "\n")
        	continue

        # Call major allele and get intervening sequence
        major_allele = vcf.get_major_allele(vcf_line, vcf)    
        seq =  vcf.sequence_between_SNPS(fa, vcf_line['CHROM'], previous_base['POS'], vcf_line['POS'])
        
        # Write to Fasta
        if previous_base['POS'] != None:
        	sequence += major_allele.lower()
        sequence += seq

        if len(sequence) / 250.0 > 1.0:
        	parts = split(sequence, 100)
        	[modified_fasta.write(p + "\n") for p in parts[:-1]]
        	sequence = parts[-1]

        # Update previous base information
        previous_base['CHROM'] = vcf_line['CHROM']
    	previous_base['POS'] = vcf_line['POS']
    	previous_base['SEQ'] = major_allele

    if line_count == None: continue
    if line_count > 5: break

    line_count += 1

modified_fasta.close()
