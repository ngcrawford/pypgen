fin = open("/Users/ngcrawford/Desktop/anoMar/fstats/Dest.fstats.csv",'rU')
fout = open("/Users/ngcrawford/Desktop/anoMar/circos/data/Dest.no_filter.fstats.txt",'w')
for count, line in enumerate(fin):
	if count == 0: continue # skip header
	chrm, start, depth, value = line.strip().split(',')
	# if float(value) < 0.2: continue
	if float(depth) < 15: continue
	if value == 'nan': continue
	stop = int(start) + 1
	final_line = "{0} {1} {2} {3}\n".format(chrm, start, stop, value)
	fout.write(final_line)