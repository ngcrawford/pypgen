
import os
import re
import sys
from collections import OrderedDict

def make_chunks(l, n):
    """ Yield n successive chunks from l.
    """
    newn = int(1.0 * len(l) / n + 0.5)
    for i in xrange(0, n-1):
        yield l[i*newn:i*newn+newn]
    yield l[n*newn-newn:]

def get_header_sizes(fin):

    chrms_sizes_dict = OrderedDict()

    for line in fin:

        if line.startswith("##contig"):
            chrm_name = re.findall(r'ID=.*,', line)
            chrm_name = chrm_name[0].strip('ID=').strip(',')

            chrm_length = re.findall(r'length=.*>', line)
            chrm_length = int(chrm_length[0].strip('length=').strip('>'))

            chrms_sizes_dict[chrm_name] = chrm_length

    fin.close()
    return chrms_sizes_dict


window = 5000
fin = open("test_header_lines.txt",'rU')
chrm_info = get_header_sizes(fin)

slices = []
count = 0
for count, i in enumerate(chrm_info.iteritems()):

    key, value = i
    print key, value
    if value < window:
        continue


    starts = range(1, value, window)[:-1]
    stops = range(0, value, window)[1:]
    chrms = [key] * len(stops)

    slices.extend(zip(chrms, starts, stops))

n = 20
l = slices
chunks = make_chunks(l, n)

for count, c in enumerate(chunks):
    fout = open("heli_5kb_regions_{}.txt".format(count),'w')
    for r in c:
        fout.write("{0}:{1}-{2}\n".format(*r))





