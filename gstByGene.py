
import pandas
import numpy as np
fin = open("CAP_MAR.Gsts.csv",'rU')

results = []


for count, line in enumerate(fin):
	line_parts = line.strip().split(",")
	chrm, pos, refbase = line_parts[:3] 

	min_coverage, max_coverage, Hs_est, Ht_est, \
	Gst_est, G_prime_st_est, G_double_prime_st_est, \
	D_est = [float(value) for value in line_parts[3:]]

	if min_coverage < 5.0: continue
	if max_coverage > 50: continue

	results.append(line_parts)
	# if count > 1000: break

columns=["chrm", "pos", "refbase", "min_coverage", 
			"max_coverage", "Hs_est", "Ht_est", "Gst_est", 
			"G_prime_st_est", "G_double_prime_st_est", "D_est"]
results = pandas.DataFrame(results, columns=columns)

s = ["min_coverage", 
			"max_coverage", "Hs_est", "Ht_est", "Gst_est", 
			"G_prime_st_est", "G_double_prime_st_est", "D_est"]

r = results.set_index(["chrm", "pos"])
r = r[s]
print r
r = r.astype(float)
print r.mean()
print r.std()


fin.close()
