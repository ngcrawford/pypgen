import VCF
import unittest
from collections import OrderedDict

class TestSlicing(unittest.TestCase):

    def setUp(self):
        self.bgzip_path = "test_data/butterfly.vcf.gz"

    def test_make_slices_default_settings(self):
        """Test slicing function with default settings: 500 bp slices"""
 
        for count, i in enumerate(VCF.get_slice_indicies(self.bgzip_path, regions=None, window_size=500)):
            if count > 10: break
        
        self.assertEqual(i, ('Chr01', 5501, 6000))
    
    def test_make_slices_default_with_params_set(self):
        """Test slicing function with window_size set"""

        for count, i in enumerate(VCF.get_slice_indicies(self.bgzip_path, regions=None, window_size=1008)):
            if count > 10: break
        
        self.assertEqual(i, ('Chr01', 11089, 12096))

    # def test_make_vcf_slices(self):
    #     print VCF.slice_vcf(self.bgzip_path, 'Chr01', 5501, 6000)

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

class TestGenoTypeParsing(unittest.TestCase):

    def test_genotypes(self):
        homo_ref = VCF.process_snp_call('0/0:10,9:19:99:254,0,337', 'A', 'T', IUPAC_ambiguities=True)
        self.assertEqual(homo_ref, 'A')

        heterozygote = VCF.process_snp_call('0/1:10,9:19:99:254,0,337', 'A', 'T', IUPAC_ambiguities=True)
        self.assertEqual(heterozygote, 'W')

        homo_alt =  VCF.process_snp_call('1/1:10,9:19:99:254,0,337', 'A', 'T', IUPAC_ambiguities=True)
        self.assertEqual(homo_alt, 'T')

        second_alt = VCF.process_snp_call('0/2:10,9:19:99:254,0,337', 'A', 'T,G', IUPAC_ambiguities=True)
        self.assertEqual(second_alt, 'R')

        double_alt = VCF.process_snp_call('1/2:10,9:19:99:254,0,337', 'A', 'T,G', IUPAC_ambiguities=True)
        self.assertEqual(double_alt ,'K')

        


if __name__ == '__main__':
    unittest.main()