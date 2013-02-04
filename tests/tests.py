#!/usr/bin/env python
# encoding: utf-8

import os
import sys
sys.path.insert(0, os.path.abspath('..'))  # Seriously?! This is fucking ugly.


import unittest
from pypgen.parser import VCF
from pypgen.misc.helpers import *
from collections import OrderedDict


class TestSlicing(unittest.TestCase):

    def setUp(self):
        self.bgzip_path = "pypgen/data/example.vcf.gz"

    def test_make_slices_default_settings(self):
        """Test slicing function with default settings: 500 bp slices"""

        for count, i in enumerate(VCF.get_slice_indicies(self.bgzip_path, regions=None, window_size=500)):
            if count > 10:
                break

        self.assertEqual(i, ('Chr01', 5501, 6000))

    def test_make_slices_default_with_params_set(self):
        """Test slicing function with window_size set"""

        for count, i in enumerate(VCF.get_slice_indicies(self.bgzip_path, regions=None, window_size=1008)):
            if count > 10:
                break

        self.assertEqual(i, ('Chr01', 11089, 12096))

    # def test_make_vcf_slices(self):
    #     print VCF.slice_vcf(self.bgzip_path, 'Chr01', 5501, 6000)

    # def test_check_length_of_first_ten_slices(self):
    #     pass

class TestVCFInfoParsing(unittest.TestCase):

    def setUp(self):
        self.bgzip_path = "pypgen/data/example.vcf.gz"

        self.populations_list = [
            "cydno:c511,c512,c513,c514,c515,c563,c614,c630,c639,c640",
            "outgroups:h665,i02-210",
            "melpo:m523,m524,m525,m589,m675,m676,m682,m683,m687,m689",
            "pachi:p516,p517,p518,p519,p520,p591,p596,p690,p694,p696"]

        self.header_dict = OrderedDict([('CHROM', None), ('POS', None),\
            ('ID', None), ('REF', None), ('ALT', None), ('QUAL', None), \
            ('FILTER', None), ('INFO', None), ('FORMAT', None), ('c511', None), \
            ('c512', None), ('c513', None), ('c514', None), ('c515', None), \
            ('c563', None), ('c614', None), ('c630', None), ('c639', None), \
            ('c640', None), ('h665', None), ('i02-210', None), ('m523', None), \
            ('m524', None), ('m525', None), ('m589', None), ('m675', None), \
            ('m676', None), ('m682', None), ('m683', None), ('m687', None), \
            ('m689', None), ('p516', None), ('p517', None), ('p518', None), \
            ('p519', None), ('p520', None), ('p591', None), ('p596', None), \
            ('p690', None), ('p694', None), ('p696', None)])

        self.allle_counts = {'melpo': {0: 12.0, 1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0},
                            'cydno': {0: 20.0, 1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0},
                            'pachi': {0: 18.0, 1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0},
                            'outgroups': {0: 2.0, 1: 2.0, 2: 0.0, 3: 0.0, 4: 0.0}}

    def test_population_string_parsing(self):
        populations = VCF.parse_populations_list(self.populations_list)

        self.assertEqual(populations, {'melpo': ['m523', 'm524', 'm525',
            'm589', 'm675', 'm676', 'm682', 'm683', 'm687', 'm689'],
            'pachi': ['p516', 'p517', 'p518', 'p519', 'p520', 'p591',
            'p596', 'p690', 'p694', 'p696'], 'cydno': ['c511',
            'c512', 'c513', 'c514', 'c515', 'c563', 'c614', 'c630',
            'c639', 'c640'], 'outgroups': ['h665', 'i02-210']})

    def test_header_to_ordered_dict_parsing(self):
        header = VCF.make_empty_vcf_ordered_dict(self.bgzip_path)

        self.assertEqual(header, self.header_dict)

    def test_header_vs_population_sample_ids(self):
        """Check that the sample IDs parsed from the population arguement
            match those in the VCF file.

            NOTE: In practice the populations arguement can contain fewer
            samples and populations than actually contained in the VCF file.
        """

        header = VCF.make_empty_vcf_ordered_dict(self.bgzip_path, )
        header_sample_ids = [item for count, item in enumerate(header) if count >= 9]

        populations_dict = VCF.parse_populations_list(self.populations_list)
        populations_sample_ids = [i for l in populations_dict.values() for i in l]

        # Check both unique IDs and equal length
        self.assertEqual(set(header_sample_ids), set(populations_sample_ids))
        self.assertEqual(len(header_sample_ids), len(populations_sample_ids))


class TestGenoTypeParsing(unittest.TestCase):

    def setUp(self):
        self.vcf_line = OrderedDict([
                        ('CHROM', 'Chr01'),
                        ('POS', 4984),
                        ('ID', '.'),
                        ('REF', 'C'),
                        ('ALT', 'T'),
                        ('QUAL', 235.91),
                        ('FILTER', 'PASS'),
                        ('INFO', {'AC': '1'}),
                        ('FORMAT', 'GT:AD:DP:GQ:PL'),
                        ('out1', {'GT': '0/0', 'GQ': 51.17, 'AD': (17, 0), 'DP': 17, 'PL': (0, 51, 669)}),
                        ('out2', {'GT': '0/0', 'GQ': 45.15, 'AD': (16, 0), 'DP': 16, 'PL': (0, 45, 600)})])

        self.populations = {'Outgroup': ['out1', 'out2']}

    #### WITH AMBIGUITY CODES ####
    def test_homo_ref_genotype_calling(self):
        homo_ref = VCF.process_snp_call('0/0:10,9:19:99:254,0,337', 'A', 'T', IUPAC_ambiguities=True)
        self.assertEqual(homo_ref, 'A')

    def test_heterozygote_genotype_calling(self):
        heterozygote = VCF.process_snp_call('0/1:10,9:19:99:254,0,337', 'A', 'T', IUPAC_ambiguities=True)
        self.assertEqual(heterozygote, 'W')

    def test_homo_alt_genotype_calling(self):
        homo_alt = VCF.process_snp_call('1/1:10,9:19:99:254,0,337', 'A', 'T', IUPAC_ambiguities=True)
        self.assertEqual(homo_alt, 'T')

    def test_second_alt_genotype_calling(self):
        second_alt = VCF.process_snp_call('0/2:10,9:19:99:254,0,337', 'A', 'T,G', IUPAC_ambiguities=True)
        self.assertEqual(second_alt, 'R')

    def test_double_alt_genotype_calling(self):
        double_alt = VCF.process_snp_call('1/2:10,9:19:99:254,0,337', 'A', 'T,G', IUPAC_ambiguities=True)
        self.assertEqual(double_alt, 'K')

    #### WITHOUT AMBIGUITY CODES ####
    def test_heterozygote_as_N_genotype_calling(self):
        heterozygote_as_N = VCF.process_snp_call('0/1:10,9:19:99:254,0,337', 'A', 'T', IUPAC_ambiguities=False)
        self.assertEqual(heterozygote_as_N, 'N')

    def test_double_alt_het_as_N_genotype_calling(self):
        double_alt_het_as_N = VCF.process_snp_call('1/2:10,9:19:99:254,0,337', 'A', 'T,G', IUPAC_ambiguities=False)
        self.assertEqual(double_alt_het_as_N, 'N')

    #### TEST OUTGROUP CALLING ####
    def test_home_ref_outgroup_calling(self):
        homo_ref = VCF.process_outgroup(self.vcf_line, self.populations)
        self.assertEqual(homo_ref, '0')

    def test_home_alt_outgroup_calling(self):
        vcf_line = self.vcf_line

        vcf_line['out1']['GT'] = '1/1'
        vcf_line['out2']['GT'] = '1/1'

        homo_alt = VCF.process_outgroup(vcf_line, self.populations)
        self.assertEqual(homo_alt, '1')

    def test_heterozyogote_outgroup_calling(self):
        vcf_line = self.vcf_line

        vcf_line['out1']['GT'] = '1/1'
        vcf_line['out2']['GT'] = '1/0'

        het = VCF.process_outgroup(vcf_line, self.populations)
        self.assertEqual(het, None)

    def test_diff_homozygotes_outgroups_calling(self):
        vcf_line = self.vcf_line

        vcf_line['out1']['GT'] = '1/1'
        vcf_line['out2']['GT'] = '0/0'

        diff = VCF.process_outgroup(vcf_line, self.populations)
        self.assertEqual(diff, None)


class TestFstatsCalculations(unittest.TestCase):

    def setUp(self):
        self.trivial_allele_counts = {'pop1': {0: 2.0, 1: 8.0, },
                                      'pop2': {0: 8.0, 1: 2.0, }}

        self.normal_allele_counts = {'melpo': {0: 18.0, 1: 2.0, 2: 0.0, 3: 0.0, 4: 0.0},
                                     'cydno': {0: 16.0, 1: 4.0, 2: 0.0, 3: 0.0, 4: 0.0},
                                     'pachi': {0: 4.0, 1: 16.0, 2: 0.0, 3: 0.0, 4: 0.0},
                                     'outgroups': {0: 4.0, 1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0}}

    def test_calc_fstats_trivial_allele_counts(self):
        f_statistics = VCF.calc_fstats(self.trivial_allele_counts)
        self.assertEqual(
            {('pop2', 'pop1'):
                {'G_prime_st_est': 0.6803049722304385,
                'D_est': 0.517460317460318,
                'G_double_prime_st_est': 0.7609710550887027,
                'Gst_est': 0.33747412008281613,
                'Hs_est': 0.3368421052631577,
                'Ht_est': 0.508421052631579}},
            f_statistics)

    def test_calc_fstats_normal_allele_counts(self):
        f_statistics = VCF.calc_fstats(self.normal_allele_counts)
        self.assertEqual(
            {'G_prime_st_est': 0.9140159767610748,
            'D_est': 0.7581699346405228,
            'G_double_prime_st_est': 0.9477124183006534,
            'Gst_est': 0.6444444444444445,
            'Hs_est': 0.1729729729729729,
            'Ht_est': 0.48648648648648646},
            f_statistics[('pachi', 'outgroups')])


class TestMultilocusFstatsCalculations(unittest.TestCase):
    def setUp(self):

        self.trivial_Hs_est_dict = {('pop1', 'pop2'): [0.1, 0.1, ]}
        self.trivial_Ht_est_dict = {('pop1', 'pop2'): [0.1, 0.1, ]}

        self.trivial_2_Hs_est_dict = {('pop1', 'pop2'): [0.3, 0.2, ]}
        self.trivial_2_Ht_est_dict = {('pop1', 'pop2'): [0.5, 0.5, ]}

        self.Hs_est_dict = {('pop1', 'pop2'): [0.0, 0.1403846153846153, ]}
        self.Ht_est_dict = {('pop1', 'pop2'): [0.14891304347826081, 0.052536231884058045, ]}

    def test_calc_multilocus_f_statistics(self):
        ml_stats = VCF.calc_multilocus_f_statistics(self.Hs_est_dict, self.Ht_est_dict)
        self.assertEqual({('pop1', 'pop2'):
            {'Gst_est': 0.30312672938572266,
            'Gst_est.stdev': 1.3360742705570265,
            'G_double_prime_st_est': 0.5003506191636901,
            'G_double_prime_st_est.stdev': 2.394045897494315,
            'G_prime_st_est': 0.34889353651117816,
            'G_prime_st_est.stdev': 1.6091544573608876,
            'D_est': 0.0015629851350084287,
            'D_est.stdev': 0.25110803099568757}}, ml_stats)

    def test_trivial_calc_multilocus_f_statistics(self):
        ml_stats = VCF.calc_multilocus_f_statistics(self.trivial_Hs_est_dict, self.trivial_Ht_est_dict)
        self.assertEqual({('pop1', 'pop2'):
            {'G_prime_st_est': 0.0,
            'Gst_est.stdev': 0.0,
            'G_double_prime_st_est.stdev': 0.0,
            'G_double_prime_st_est': 0.0,
            'Gst_est': 0.0,
            'D_est.stdev': 0.0,
            'G_prime_st_est.stdev': 0.0,
            'D_est': 0.0}}, ml_stats)

    def test_trivial_2_calc_multilocus_f_statistics(self):
        ml_stats = VCF.calc_multilocus_f_statistics(self.trivial_2_Hs_est_dict, self.trivial_2_Ht_est_dict)
        self.assertEqual({('pop1', 'pop2'):
        {'G_prime_st_est': 0.8333333333333334,
        'Gst_est.stdev': 0.09999999999999998,
        'G_double_prime_st_est.stdev': 0.060586734693877375,
        'G_double_prime_st_est': 0.8888888888888888,
        'Gst_est': 0.5,
        'D_est.stdev': 0.08928571428571419,
        'G_prime_st_est.stdev': 0.07857142857142851,
        'D_est': 0.6488650338184054}}, ml_stats)


if __name__ == '__main__':
    unittest.main()
