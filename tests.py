import VCF
import unittest
from collections import OrderedDict

class TestSlicing(unittest.TestCase):

    def setUp(self):
        self.bgzip_path = "test_data/butterfly.vcf.gz"

    def test_make_slices_default_settings(self):
        """Test slicing function with default settings: 500 bp slices"""
 
        for count, i in enumerate(VCF.get_slice_indicies(self.bgzip_path)):
            if count > 10: break
        self.assertEqual(i, ('Chr01', 5501, 6000))
    
    def test_make_slices_default_with_params_set(self):
        """Test slicing function with window_size set"""

        for count, i in enumerate(VCF.get_slice_indicies(self.bgzip_path, window_size=1008)):
            if count > 10: break
        self.assertEqual(i, ('Chr01', 11089, 12096))

    def test_make_vcf_slices(self):
        print VCF.slice_vcf(self.bgzip_path, 'Chr01', 5501, 6000)

    def test_check_length_of_first_ten_slices(self):
        pass

class TestVCFInfoParsing(unittest.TestCase):

    def setUp(self):
        self.bgzip_path = "test_data/butterfly.vcf.gz"
        
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

        self.vcf_line = OrderedDict([('CHROM', 'Chr01'), 
                                     ('POS', '5906'), 
                                     ('ID', '.'), 
                                     ('REF', 'G'), 
                                     ('ALT', 'A'), 
                                     ('QUAL', '7211.87'), 
                                     ('FILTER', 'PASS'), 
                                     ('INFO', {'AC': '22', 
                                               'BaseQRankSum': '0.013', 
                                               'MQRankSum': '-3.442', 
                                               'InbreedingCoeff': '0.4448', 
                                               'AF': '0.344', 
                                               'HRun': '2', 
                                               'AN': '64', 
                                               'MQ0': '1', 
                                               'Dels': '0.01', 
                                               'FS': '11.310', 
                                               'MQ': '63.35', 
                                               'QD': '27.21', 
                                               'HaplotypeScore': '3.6614', 
                                               'SB': '-2692.78', 
                                               'DP': '535', 
                                               'ReadPosRankSum': 
                                               '-0.922'}), 
                                     ('FORMAT', 'GT:AD:DP:GQ:PL'), 
                                     ('c511', {'GT': '0/0', 'GQ': '48.15', 'AD': '19,0', 'DP': '19', 'PL': '0,48,622'}), 
                                     ('c512', {'GT': '0/0', 'GQ': '45.11', 'AD': '15,0', 'DP': '15', 'PL': '0,45,553'}), 
                                     ('c513', {'GT': '0/0', 'GQ': '69.20', 'AD': '25,0', 'DP': '25', 'PL': '0,69,886'}), 
                                     ('c514', {'GT': '0/0', 'GQ': '57.18', 'AD': '21,0', 'DP': '21', 'PL': '0,57,727'}), 
                                     ('c515', {'GT': '0/1', 'GQ': '63.50', 'AD': '3,8', 'DP': '11', 'PL': '245,0,64'}), 
                                     ('c563', {'GT': '0/1', 'GQ': '99', 'AD': '8,11', 'DP': '19', 'PL': '327,0,243'}), 
                                     ('c614', {'GT': '0/1', 'GQ': '99', 'AD': '8,12', 'DP': '20', 'PL': '383,0,253'}), 
                                     ('c630', {'GT': '0/0', 'GQ': '57.18', 'AD': '20,0', 'DP': '20', 'PL': '0,57,745'}), 
                                     ('c639', {'GT': '0/1', 'GQ': '99', 'AD': '9,9', 'DP': '18', 'PL': '277,0,109'})
  
    def test_population_string_parsing(self):
        populations = VCF.parse_populations_list(self.populations_list)
        self.assertEqual(populations, {'melpo': ['m523', 'm524', 'm525', \
            'm589', 'm675', 'm676', 'm682', 'm683', 'm687', 'm689'], \
            'pachi': ['p516', 'p517', 'p518', 'p519', 'p520', 'p591', \
            'p596', 'p690', 'p694', 'p696'], 'cydno': ['c511', \
            'c512', 'c513', 'c514', 'c515', 'c563', 'c614', 'c630', \
            'c639', 'c640'], 'outgroups': ['h665', 'i02-210']})

    def test_header_to_ordereddict_parsing(self):
        header = VCF.set_header(self.bgzip_path)
        self.assertEqual(header,self.header_dict)

    def test_header_vs_population_sample_ids(self):
        """Check that the sample IDs parsed from the population arguement
            match those in the VCF file.

            NOTE: In practice the populations arguement can contain fewer 
            samples and populations than actually contained in the VCF file.
        """

        header = VCF.set_header(self.bgzip_path)
        header_sample_ids = [item for count, item in enumerate(header) if count >= 9]
        
        populations_dict  = VCF.parse_populations_list(self.populations_list)
        populations_sample_ids = [i for l in populations_dict.values() for i in l]

        # Check both unique IDs and equal length
        self.assertEqual(set(header_sample_ids), set(populations_sample_ids))
        self.assertEqual(len(header_sample_ids), len(populations_sample_ids))


if __name__ == '__main__':
    unittest.main()