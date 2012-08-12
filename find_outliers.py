#!/usr/bin/env python
# encoding: utf-8

"""
find_outliers.py

This is just an initial pass at a method/script. 
"""

import os
import glob
import pandas


def find_outliers(fname, proportion=0.05):

	fout = os.path.splitext(fname)[0] + '{0}.circos.txt'.format(proportion)
	d = pandas.read_csv(fname,sep="\t")
	dn= d.dropna()
	stop = int(dn.shape[0] * proportion)
	dns = dn.sort(columns='Gst_est', axis=0, ascending=False)
	dns = dns[["CHROM","POS",'Gst_est']]
	print stop, dns[['Gst_est']].values[stop]
	stop = dns[['POS']].values.T[0] + 1
	dns.insert(2, 'STOP', stop)
	dns = dns.rename(columns={'POS':'START'})
	dns.to_csv(fout, sep=' ', header=False, index=False)

def ht_vs_Gst(fname):
	d = pandas.read_csv(fname,sep="\t")
	z = d.dropna()

	Ht = z[['Ht_est']].values.T[0]
	Gst_est = z[['Gst_est']].values.T[0]

	xlabel('Ht_est')
	ylabel('Gst_est')
	title('Marina vs Plage de Clugny')
	figsize(10,10)

def remove_zeros(fname):
	fin = open(fname,'rU')
	fout = open("/Users/ngcrawford/Desktop/anoMar/fstats/MAR_CAP.fstats.stampy.5_samples_min.noZeros.txt", 'w')
	for count, line in enumerate(fin):
		if count % 100000 == 0:
			print count, line
		line_parts = line.strip().split()
		if line_parts[-1] == '0.0':
			continue
		else:
			fout.write(line)


def make


def make_circos_input(fname):
	fin = open(fname,'rU')
	fout = open("/Users/ngcrawford/Desktop/anoMar/fstats/MAR_CAP.Gst_est.circos.txt", 'w')
	for count, line in enumerate(fin):
		line_parts = line.strip().split()

		if count == 0:
			continue
		
		if float(line_parts[-2]) != 0.39 or line_parts[-1] == 'nan':
			continue

		else:
			chrm = line_parts[0]
			start = line_parts[1]
			stop = str(int(start) +1)
			gst = line_parts[6]
			line = [chrm, start, stop, gst]
			line = ' '.join(line) + '\n'
			fout.write(line)


if __name__ == '__main__':

	#for fname in glob.glob("Users/ngcrawford/Desktop/anoMar/fstats/MAR_CAP.fstats.stampy.5_samples_min.txt"):
	make_circos_input("/Users/ngcrawford/Desktop/anoMar/fstats/MAR_CAP.fstats.stampy.5_samples_min.noZeros.txt")



