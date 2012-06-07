fin = open("/Users/ngcrawford/Desktop/anoMar/circos/data/Amar.CAP-MAR.Gst_ext.w1000bp.csv",'rU')
fout = open("/Users/ngcrawford/Desktop/anoMar/circos/data/Amar.CAP-MAR.Gst_est.w1000bp.filtered.csv",'w')
for count, line in enumerate(fin):
	line_parts = line.strip().split(',')
	if len(line_parts) < 5: continue
	if count == 0: continue # skip header
	chrm, start, stop, depth, value = line.strip().split(',')
	if float(value) < 0.2: continue
	if value == 'nan': continue
	final_line = "{0} {1} {2} {3}\n".format(chrm, start, stop, value)
	fout.write(final_line)