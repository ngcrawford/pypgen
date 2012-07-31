import pylab
import numpy
import pandas
import random

def find_non_overlapping_slices(data,expected_overlap = 5000):
	previous_chrm = None
	previous_start = None
	non_overlapping_data = []
	row_count = 0
	for count, row in enumerate(data.iterrows()):
	    #if count > 100: break
	    if previous_chrm == None:
	        previous_chrm = row[1]['chrm']
	        previous_start = row[1]['start']
	    
	    if previous_chrm != row[1]['chrm']:
	        previous_chrm = row[1]['chrm']
	        previous_start = row[1]['start']
	        continue
	    
	    if row[1]['start'] - previous_start >= expected_overlap:
	        non_overlapping_data.append(row[1].to_dict())
	        previous_chrm = row[1]['chrm']
	        previous_start = row[1]['start']
	        row_count += 1

	return pandas.DataFrame(non_overlapping_data)

def calc_pvalues(measured_values, simulated_values):
    """ Given a list of measured values and a list of simulated values tailed probabilites (e.g., p-values) are calculated"""
    
    if type(simulated_values) == pandas.core.series.Series:
    	simulated_values = simulated_values.values

    if type(measured_values) == pandas.core.series.Series:
    	measured_values = measured_values.values

    simulated_values.sort()
    simulated_values = simulated_values[~numpy.isnan(simulated_values)]
    inserts = numpy.searchsorted(simulated_values, measured_values)
    pvalues = numpy.array([(len(simulated_values) - float(item)) / len(simulated_values) for item in inserts])
    return pvalues

def generate_random_sample(data, draws=10000, cols_2_skip = ['chrm', 'start', 'stop']):

	results = {}
	columns = [col for col in non_overlapping_data.columns if col not in cols_2_skip]
	for count, col in enumerate(columns):
		rand_values = numpy.array([random.choice(data[col].values) for item in xrange(0,draws)])
		results[col] = rand_values

	rand_sample = pandas.DataFrame(results)
	return rand_sample

def adjust_pvalues(pvalues):
	from rpy2.robjects.packages import importr
	from rpy2.robjects.vectors import FloatVector
	stats = importr('stats')
	p_adjust = numpy.array(stats.p_adjust(FloatVector(pvalues), method = 'BH',))
	return p_adjust

def update_data_with_pvalues(data, simulated_data, cols_2_skip = ['chrm', 'start', 'stop']):
	columns = [col for col in non_overlapping_data.columns if col not in cols_2_skip]
	for col in columns:
		name = col + ".pvalues"
		pv = calc_pvalues(data[col], simulated_data[col])
		data[name] = pv

	return data





fin = open('/Users/ngcrawford/Desktop/anoMar/circos/data/sliding_window_stampy/CAP_MAR.5000bp_w.2500bp_ol.dadi.corrected_headers.txt','rU')
data = pandas.read_csv(fin)
print 'read data'
non_overlapping_data = find_non_overlapping_slices(data)
print 'made non_overlapping_data'
random_sample = generate_random_sample(non_overlapping_data)
print 'made random sample'
data = update_data_with_pvalues(data, random_sample)

hist(data['Fst.pvalues'])

