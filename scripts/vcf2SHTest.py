#!/usr/bin/env python
# encoding: utf-8

"""
vcf2trees.py

Created by Nick Crawford on 2011-11-18.
Copyright (c) 2011

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses

The author may be contacted at ngcrawford@gmail.com
"""

import os
import re
import sys
import gzip
import glob
import pysam
import shlex
import pandas
import random
import unittest
import tempfile
import argparse
import itertools
import numpy as np
import pypgen
from pypgen.parser import VCF
from copy import deepcopy, copy
from subprocess import Popen, PIPE
from collections import namedtuple, defaultdict

"""
python scripts/vcf2SHTest.py \
-r Chr01:0-5000 \
pypgen/data/example.vcf.gz \
-c "((((c511,c512,c513,c514,c515,c563,c614,c630,c639,c640),(p516,p517,p518,p519,p520,p591,p596,p690,p694,p696)),(m523,m524,m525,m589,m675,m676,m682,m683,m687,m689)),(h665,i02_210));((((c511,c512,c513,c514,c515,c563,c614,c630,c639,c640),(m523,m524,m525,m589,m675,m676,m682,m683,m687,m689)),(p516,p517,p518,p519,p520,p591,p596,p690,p694,p696)),(h665,i02_210));((((m523,m524,m525,m589,m675,m676,m682,m683,m687,m689),(p516,p517,p518,p519,p520,p591,p596,p690,p694,p696)),(c511,c512,c513,c514,c515,c563,c614,c630,c639,c640)),(h665,i02_210));" \
--sh-test \
-p 3
"""

def get_args():
    """Parse sys.argv"""
    parser = argparse.ArgumentParser()

    parser.add_argument('-r', '-R', '--regions',
                        required=True,
                        type=str,
                        help="Chromosomal region in the format: 'chrm:start-stop'")

    parser.add_argument('-o', '--output',
                        type=argparse.FileType('w'),
                        default=sys.stdout,
                        help='Path to output. (default is STOUT)')

    parser.add_argument('-c', '--constraint-trees',
                        type=str,
                        help="Newick formated trees.")

    parser.add_argument('-m', '--model',
                        type=str,
                        default='GTRGAMMA',
                        help="Substitution model name. \
                              GTRGAMMA (default) | GTRCAT | GTRCATI | GTRGAMMAI")

    # parser.add_argument('-b','--bootstraps',
    #                     type=int,
    #                     help='Calculate bootstraps.')

    parser.add_argument('-p', '--threads',
                        default=1,
                        type=int,
                        help='Set number of threads to use. Only for RAXML.')

    parser.add_argument('--raxml',
                        action='store_true',
                        default=False,
                        help="beta: run raxml in place of phyml. Works with constraint trees.")

    # parser.add_argument('--name',
    #                     type=str,
    #                     help='Add uniqe id to output.')

    parser.add_argument('--as-nexus',
                        action="store_true",
                        default=False)

    parser.add_argument('--sh-test',
                        action="store_true",
                        default=False)

    #parser.add_argument('--as-alignment')

    parser.add_argument('input',
                        nargs=1,
                        help='Bgzipped and indexed VCF file')

    args = parser.parse_args()
    return args

class SHTest(object):
    """docstring for SHTest"""
    def __init__(self):
        super(SHTest, self).__init__()

    def count_informative_sites(self, alignment_array):
        """Informative Sites must have two different SNPs"""
        informative_sites = 0

        for site in alignment_array:
            unique_sites = set(site)

            if len(unique_sites) >= 3:
                informative_sites += 1

        return informative_sites

    def array2OnelinerAlignment(self, info, taxa, bases):
        """Convert array of array of taxa and an array of bases to one-liner."""

        oneliner = info

        for count, seq in enumerate(bases):
            oneliner += taxa[count] + "," + ''.join(itertools.chain(bases[count])) + ","

        oneliner = oneliner[:-1] + ";"
        return oneliner

    def makeDataTuple(self, vcf):
        """Setup a labeled tuple to store the data."""

        chrm_data = {}

        for count, line in enumerate(vcf):

            if "##contig" in line:
                contigline = line.split("<")[1]
                contigline = contigline.strip(">\n").split(",")
                contigline = [item.split("=")[1] for item in contigline]
                chrm_data[contigline[0]] = int(contigline[1])

            if line.startswith("#CHROM"):
                field_labels = line.strip().strip("#").split("\t")
                field_labels = [item.strip().split('.')[0].replace("-", "_") for item in field_labels]
                break

        position_data = namedtuple('base', field_labels)

        return (position_data, chrm_data)

    def process_vcf_slice(self, tabix_file, region, position_data):

        chrm, start, stop = region['chrm'], region['start'], region['stop']

        tbx = pysam.Tabixfile(tabix_file)
        tbx_lines = tbx.fetch(chrm, start, stop)

        numb_of_seqs = len(position_data._fields[9:])
        alignment = np.zeros((stop - start, numb_of_seqs), np.string0)

        # This 'error handling' needs to be rewritten.
        current_data = []

        if tbx_lines is None:
            return 'error'

        for line in tbx_lines:
            current_base = position_data._make(line.strip().split("\t"))
            base_calls = self.callSNPs(current_base, numb_of_seqs)
            current_data.append(base_calls.copy())

        alignment = np.array(current_data)
        inform_sites = self.count_informative_sites(alignment)

        if current_base is None:
            return 'error'
        else:
            taxa = current_base._fields[9:]
            info = "tree 'chrm={0},start={1},stop={2},inform_sites={3}':".format(current_base.CHROM, start, stop, inform_sites)
            oneliner = self.array2OnelinerAlignment(info, taxa, alignment.transpose())

        if ":" in oneliner and oneliner[-1] == ';':   # this prevents bad alignments from getting printed
            return oneliner
        else:
            return 'error'

    def callSNPs(self, current_base, numb_of_seqs):
        """Call the SNPs. Duh!"""

        blanks = np.zeros(numb_of_seqs, np.string0)

        if current_base.FILTER == 'LowQual':
            blanks.fill("-")

        if current_base.FORMAT == 'GT':
            blanks.fill("-")

        for count, snp_call in enumerate(current_base[9:]):
            base = VCF.process_snp_call(snp_call, current_base.REF, current_base.ALT)
            blanks[count] = base

        return blanks

    def oneliner2phylip(self, line):
        """Convert one-liner to phylip format."""

        seqs = line.strip(";").split(':')[-1].split(',')
        label_seqs = zip(seqs[:-1:2], seqs[1::2])
        taxa_count = len(label_seqs)
        seq_length = len(label_seqs[0][1])
        alignment = "%s %s\n" % (taxa_count, seq_length)   # add header

        for taxa_name, seq in label_seqs:
            taxa_name = taxa_name.strip()
            alignment += '%-10s%s\n' % (taxa_name, seq)

        return alignment

    def region_2_dict(self, args):
        """ Convert region string to dictionary.
            e.g., Chr01:0-5000 => {'start': 0, 'chrm': 'Chr01', 'stop': 5000} """

        chrm, start, stop = re.split(r':|-', args.regions)
        start, stop = int(start), int(stop)
        return {'chrm': chrm, 'start': start, 'stop': stop}

    def calculate_raxml_trees(self, phylip, model, threads, pos, const_tree=None):
        """Doc string"""

        args_dict = pos
        if os.path.exists('tmp/') is False:
            os.mkdir('tmp/')

        temp_in = tempfile.NamedTemporaryFile(suffix='.out', dir='tmp/')   # delete=False)
        for line in phylip:
            temp_in.write(line)
        temp_in.seek(0)     # move pointer to beginning of file

        temp_string = os.path.split(temp_in.name)[1].split('.')[0]

        # SETUP CONSTRAINT TREE FILE
        if const_tree is not None:

            constraint_file = tempfile.NamedTemporaryFile(suffix='.tree', dir='tmp/')
            constraint_file.write(const_tree + '\n')
            constraint_file.seek(0)

            # raxmlHPC-SSE3 -m GTRGAMMA -s test.phylip -n t -g test.constraint.tree -o h665,i02_210

            cli = 'raxmlHPC-PTHREADS-SSE3  \
                    -m {1} \
                    -s {0} \
                    -g {2} \
                    -n {3} \
                    -T {4} \
                    -o h665,i02_210 \
                    >/dev/null 2>&1'.format(temp_in.name, model, constraint_file.name, temp_string, threads)

        else:
            cli = 'raxmlHPC-PTHREADS-SSE3 \
                   -m {1} \
                    -s {0} \
                    -n {2} \
                    -T {3} \
                    -o h665,i02_210 \
                    >/dev/null 2>&1'.format(temp_in.name, model, temp_string, threads)

        cli_parts = cli.split()
        ft = Popen(cli_parts, stdin=PIPE, stderr=PIPE, stdout=PIPE).communicate()

        temp_string = os.path.split(temp_in.name)[1].split('.')[0]

        treefile = 'RAxML_result.%s' % (temp_string)
        tree = open(treefile, 'r').readlines()[0].strip().strip("\"")

        statsfile = 'RAxML_info.%s' % (temp_string)
        lnL = self.processStatsFile(open(statsfile, 'r'))

        args_dict['lnL'] = float(lnL)
        args_dict['model'] = model

        args_dict['tree'] = "'" + str(tree) + "'"

         ## clean up:
        [os.remove(f) for f in glob.glob("RAxML_*")]

        return args_dict

    def processStatsFile(self, fin):
        lnL = None
        for line in fin:
            if 'Log-likelihood' in line:
                lnL = line.split()[-1]

            if line.startswith("Final GAMMA-based Score of best tree") is True:
                lnL = line.split()[-1]

        return lnL

    def calculate_sh_test(self, phylip, region, model, threads, best_tree, const_trees, fitted_trees):
        """Do SH test"""

        # WRITE TREES TO TEMP FILES
        best_file = tempfile.NamedTemporaryFile(suffix='.best', dir='tmp/')   # delete=False)
        best_file.write(best_tree["tree"] + "\n")
        best_file.seek(0)

        const_file = tempfile.NamedTemporaryFile(suffix='.const', dir='tmp/')   # delete=False)
        for t in fitted_trees:
            const_file.write(t["tree"] + "\n")
        const_file.seek(0)     # move pointer to beginning of file

        # WRITE ALIGNMENT TO TEMP FILE
        temp_in = tempfile.NamedTemporaryFile(suffix='.out', dir='tmp/')   # delete=False)
        for line in phylip:
            temp_in.write(line)
        temp_in.seek(0)     # move pointer to beginning of file

        # RUN SH TEST
        run_name = os.path.split(temp_in.name)[1].split('.')[0]


        cli = "raxmlHPC-SSE3 \
          -f h \
          -t {0}  \
          -z {1} \
          -s {2} \
          -m GTRGAMMA \
          -n {3}".format(best_file.name, const_file.name, temp_in.name, run_name)

        cli_parts = cli.split()
        ft = Popen(cli_parts, stdin=PIPE, stderr=PIPE, stdout=PIPE).communicate()

        # PARSE SH TEST OUTPUT
        SHTest_results = []
        best_tree_line = None
        for line in ft[0].split('\n'):

            if line.startswith('Tree:') is True:
                line = re.split('\s+', line)
                line = [w.replace(":", '') for w in line]
                SHTest_results.append(line)

            elif line.startswith('Model optimization') is True:
                line = re.split('\s+', line)
                line = [w.replace(":", '') for w in line]
                line += [best_tree["tree"]]
                best_tree_line = line

        # CLEAN UP
        [os.remove(f) for f in glob.glob("*.%s" % (run_name))]

        return (SHTest_results, best_tree_line)

    def deconvolve_sh_results(self, z):

        tnames = list(z['Tree_IDs'])
        sh_tests = z["SH_tree_name_results"] 
        sh_results = z["SH_test_results"]

        # This particular bit of logic was SHOCKINGLY difficult
        # to write.

        final_list = []
        final = None
        current_tnames = tnames
        for count, tlist in enumerate(sh_tests):

            for c, i in enumerate(current_tnames):

                if i not in tlist and i != 'best.tre' and i != 'NA':
                    current_tnames[c] = 'best.tre'
                    final[c] = 'best.tre'

                if i in current_tnames and i != 'best.tre' and i != 'NA':
                    z = tlist.index(i)

                    final[c] = sh_results[count][z]

                if c == 0:
                    final = copy(current_tnames)

            current_tnames = ['NA' if x == 'best.tre' else x for x in current_tnames]
            final_list.append(final)

        return final_list

    def recursively_sh_test(self, phylip, region, model, threads, best_tree,\
                            const_trees, fitted_trees, tree_names, results_template=None,\
                            iteration=0, best_const_tree_ids=[]):
        """
        1.) Initial test calculates best tree
            - subsequent tests use tree with largest -lnL as best tree

        2.) Repeat test until no constraint trees are left.
        """

        SHTest_results, best_tree_line = self.calculate_sh_test(phylip, region, model, threads, best_tree, const_trees, fitted_trees)

        SHTest_results = np.array(SHTest_results)
        lnLs = SHTest_results[:, 3].astype(np.float)
        best_const_tree_id = lnLs.argmax()
        best_const_tree_ids.append(best_const_tree_id)

        if iteration != 0:
            results_template['best_const_pos'] = best_const_tree_ids
            #print const_tree_names

        else:
            tree_names = tree_names[1:]

        
        results_template["SH_tree_name_results"].append(deepcopy(tree_names))
        results_template["SH_test_results"].append(deepcopy(list(SHTest_results[:, -2])))

        
        #del const_tree_names[best_const_tree_id]
        del const_trees[best_const_tree_id]
        del fitted_trees[best_const_tree_id]
        del tree_names[best_const_tree_id]

        """
        Region          treename    ?   model_evol      lnL scaffold    start   stop    tree here SH-Test SH-test of constraints
                                                                 
        HE671174:390001-395000  best.tre    TRUE    GTRGAMMA    -1800.635416    HE671174    390001  395000  TREE    NA      NA
        HE671174:390001-395000  constraint1 TRUE    GTRGAMMA    -1893.713626    HE671174    390001  395000  TREE    Yes(1%)     Yes (?)
        HE671174:390001-395000  constraint2 TRUE    GTRGAMMA    -1889.310928    HE671174    390001  395000  TREE    Yes(1%)     No (?)
        HE671174:390001-395000  constraint3 TRUE    GTRGAMMA    -1888.095646    HE671174    390001  395000  TREE    Yes(1%)     BEST_C
        """

        iteration += 1

        if len(const_trees) != 0:
            return self.recursively_sh_test(phylip, region, model, threads, best_tree, const_trees, fitted_trees, tree_names, results_template, iteration, best_const_tree_ids)
        else:
            return results_template

    def write_csv(self, z):
        del z["SH_test_results"]
        del z["best_const_pos"]
        del z["SH_tree_name_results"]
        df = pandas.DataFrame(z)
        cols = ['Region','Chrm', 'Start', 'Stop', 'Model', 'Tree_IDs', 'Trees', 'lnL', 'SH_Comparison_1', 'SH_Comparison_2', 'SH_Comparison_3']
        df = df[cols]
        fout = '{0}.csv'.format(df["Region"][0])
        df.to_csv(fout, sep="\t", header=True, index=False)

    def run_test(self):

        args = get_args()
        region = self.region_2_dict(args)

        # OPEN VCF
        vcf_file = gzip.open(args.input[0], 'rb')
        position_data, chrm_data = self.makeDataTuple(vcf_file)
        vcf_file.close()

        # CONVERT SLICE TO PHYLIP
        oneliner = self.process_vcf_slice(args.input[0], region, position_data)
        phylip = self.oneliner2phylip(oneliner)sb

        # CALCULATE BEST TREE
        best_tree = deepcopy(self.calculate_raxml_trees(phylip, args.model, args.threads, region, const_tree=None))

        # CALCUATE FITTED CONSTRAINT TREES
        fitted_trees = []
        const_trees = args.constraint_trees.strip(";").split(";")
        for count, c in enumerate(const_trees):
            c += ";"
            t = self.calculate_raxml_trees(phylip, args.model, args.threads, region, const_tree=c)
            fitted_trees.append(deepcopy(t))

        inc_count = len(const_trees) + 1
        tree_names = ['best.tre'] + ['constraint{0}'.format(i) for i in range(len(const_trees))]

        results_template = {
                            "Region": [args.regions] * inc_count,
                            "Tree_IDs": tree_names,
                            "Model": [args.model] * inc_count,
                            "lnL": [best_tree['lnL']] + [t['lnL'] for t in fitted_trees],
                            'Chrm': [region["chrm"]] * inc_count,
                            'Start': [region['start']] * inc_count,
                            'Stop': [region['stop']] * inc_count,
                            'Trees': [best_tree['tree']] + [t['tree'] for t in fitted_trees],
                            'SH_tree_name_results': [],
                            'SH_test_results': []
                            }


        updated_results = self.recursively_sh_test(phylip, region, args.model, args.threads, best_tree, const_trees, fitted_trees, tree_names, results_template)
        SH_results = self.deconvolve_sh_results(updated_results)

        for ct, comparison in enumerate(SH_results,1):
            updated_results["SH_Comparison_{0}".format(ct)] = comparison

        self.write_csv(updated_results)


        #print " ".join(['Region', "Model", '-lnl', 'Chrm', 'Start', 'Stop', 'Topology', 'Constraint', 'SHTest','SHTest_Consts'])

sh = SHTest()
sh.run_test()



# class TestSHTest(unittest.TestCase):

#     def setUp(self):
#         module_dir = os.path.dirname(pypgen.__file__)
#         self.phylip_path = os.path.join(module_dir, "data/SHTest.phylip")
#         self.phylip = ""

#         for l in open(self.phylip_path, 'rU').readlines():
#             self.phylip += l

#         self.sh = SHTest()

#     def test_raxml_calculation(self):
#         model = "GTRGAMMA"
#         threads = 3
#         region = {'chrm': 'chrm', 'start': 1, 'stop': 1000}

#         tree_and_data = self.sh.calculate_raxml_trees(self.phylip, model, 3, region, const_tree=None)
#         lnL = tree_and_data["lnL"]
#         print lnL


# unittest.main()





            
