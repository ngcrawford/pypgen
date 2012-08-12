
import gzip
import dadi
import pysam
import numpy as np
from copy import copy
from random import choice
from fstats import fstats
from itertools import combinations
from collections import OrderedDict

class VCF(object):
    """docstring for VCF"""
    def __init__(self):
        super(VCF, self).__init__()
    
        self.header = None
        self.__header_dict__ = None

        self.sample_format = None
        self.sample_format_dict = None

        self.populations = None
        self.f_statistic  = None
        self.vcf_file = None
        self.window_size = None
        self.chrm2length = None

    def parse_individual_snps(self, vcf ):

        if vcf.endswith(".gz") == True:
            vcf = gzip.open(vcf,'rb')

        else:
            vcf = open(vcf,'rU')
        

        # SETUP NAMED TUPLE TO STORE INFO FROM A SINGLE BASE
        field_labels = []

        # SETUP COUNTERS
        current_base = None
        line_count = None

        # PARSE VCF FIlE
        snp_count = 0

        for line in vcf:

            # SKIP HEADER
            if line.startswith("#CHROM"):
                self.header = line.strip("#").strip().split()
                self.__header_dict__ = OrderedDict([(item,None) for item in self.header])
                line_count = 0

            # START PROCESSING ALIGNED BASES
            if line_count > 0:
                vcf_line_dict = self.parse_vcf_line(line)
                yield vcf_line_dict

            if line_count >=0:
                line_count += 1

    def filter_vcf_line(self, filter_string, vcf_line):
        """Apply filter to VCF. If VCF passes filter returns True.

            Example filter string: "QUAL > 500"

            Should work with JEXL expressions used in GATK. But, I 
            should think about the best, most pythonic, way to do this.
        """

        # TO DO: rewrite
        # This is very unsophisticated and insecure!


        col, exp, value = filter_string.split(' ')
        exp = "vcf_line[%s] %s %s" % (col, exp, value )

        return eval(exp)


    def parse_info(self,info_field):
        
        info = []
        for item in info_field.split(','): # TO DO: comma should be ';'
            pair = item.split("=") 
            if len(pair) == 2:
                info.append(pair)

        info_dict = dict(info)
        return info_dict

    def parse_vcf_line(self, line):

        line_parts = line.strip().split()
        vcf_line_dict = self.__header_dict__

        for count, item in enumerate(vcf_line_dict):
            vcf_line_dict[item] = line_parts[count]

        if self.sample_format == None:
            self.sample_format = vcf_line_dict["FORMAT"].split(":")
            self.sample_format_dict = dict([(item,None) for item in self.sample_format])      
     
        for count, item in enumerate(vcf_line_dict):
            
            if count >= 9:
                genotype = vcf_line_dict[item]
                
                if genotype == "./.":
                    vcf_line_dict[item] = None

                else:
                    genotype = dict(zip(self.sample_format,genotype.split(":")))
                    vcf_line_dict[item] = genotype

        return vcf_line_dict

   
    def calc_fstats(self, allele_counts):

        #test_pair = {'So_Riviere_Goyaves': np.array([ 0.0, 1.0, 0.0, 0.0]), 'Plage_de_Viard': np.array([ 1.0, 0.0, 0.0, 0.0]),}
        
        # CALCULATE ALLELE FREQUENCIES
        allele_freqs_dict = self.populations.fromkeys(self.populations.keys(),None)
        for pop in self.populations:
            counts =  allele_counts[pop].values()
            freqs =  counts/np.sum(counts, dtype=float)
            allele_freqs_dict[pop] = freqs
        
        allele_freqs = allele_freqs_dict.values()


        # CACULATE PAIRWISE F-STATISTICS
        pairwise_results = {}

        pops_no_outgroup =  self.populations.keys()
        if 'outgroups' in pops_no_outgroup:
           pops_no_outgroup.remove('outgroups')

        for population_pair in combinations(pops_no_outgroup,2):
            
            pop1, pop2 =  population_pair

            Ns = [sum(allele_counts[pop].values()) for pop in [pop1, pop2]]

            if 0 in Ns: 
                values = [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]
                values_dict = dict(zip(['Hs_est', 'Ht_est', 'Gst_est', 'G_prime_st_est', 'G_double_prime_st_est', 'D_est'],values))
                pairwise_results[population_pair] = values_dict
                values_dict['CHROM'] = allele_counts["CHROM"]
                values_dict['POS'] = allele_counts["POS"]
                continue
            
            else:
 
                pop1 = allele_freqs_dict[pop1]
                pop2 = allele_freqs_dict[pop2]

                # Skip populations fixed at particular SNPS
                # this can happen when the reference is divergent.
                if pop1[1] == 1.0 and pop2[1] == 1.0:
                    values = [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]
                    values_dict = dict(zip(['Hs_est', 'Ht_est', 'Gst_est', 'G_prime_st_est', 'G_double_prime_st_est', 'D_est'],values))
                    pairwise_results[population_pair] = values_dict
                    values_dict['CHROM'] = allele_counts["CHROM"]
                    values_dict['POS'] = allele_counts["POS"]
                    continue # skip fixed SNPs
               
                if pop1[0] == 1.0 and pop2[0] == 1.0:
                    values = [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]
                    values_dict = dict(zip(['Hs_est', 'Ht_est', 'Gst_est', 'G_prime_st_est', 'G_double_prime_st_est', 'D_est'],values))
                    pairwise_results[population_pair] = values_dict
                    values_dict['CHROM'] = allele_counts["CHROM"]
                    values_dict['POS'] = allele_counts["POS"]
                    continue # skip fixed SNPs 

                print pop1, pop2


                allele_freqs = [pop1, pop2]

                n = float(len(Ns))
                Ns_harm = fstats.harmonic_mean(Ns)
                Ns_harm_chao = fstats.harmonic_mean_chao(Ns)

                # CALCULATE Hs AND Ht
                n = 2
                Hs_prime_est_ = fstats.Hs_prime_est(allele_freqs,n)
                Ht_prime_est_ = fstats.Ht_prime_est(allele_freqs,n)
                Hs_est_ = fstats.Hs_est(Hs_prime_est_, Ns_harm)
                Ht_est_ = fstats.Ht_est(Ht_prime_est_,Hs_est_,Ns_harm,n)

                # CALCULATE F-STATISTICS
                Gst_est_ = fstats.Gst_est(Ht_est_, Hs_est_)
                G_prime_st_est_ = fstats.G_prime_st_est(Ht_est_, Hs_est_, Gst_est_, n)
                G_double_prime_st_est_ = fstats.G_double_prime_st_est(Ht_est_, Hs_est_, n)
                D_est_ = fstats.D_est(Ht_est_, Hs_est_, n)
                
                # PRINT OUTPUT
                values = [Hs_est_, Ht_est_, Gst_est_, G_prime_st_est_, G_double_prime_st_est_, D_est_]
                values_dict = dict(zip(['Hs_est', 'Ht_est', 'Gst_est', \
                                        'G_prime_st_est', 'G_double_prime_st_est', 'D_est'],values))

                values_dict['CHROM'] = allele_counts["CHROM"]
                values_dict['POS'] = allele_counts["POS"]

                pairwise_results[population_pair] = values_dict

        return pairwise_results

    def __write_to_outfiles__(self, chrm, pos, depth, stat_id, f_statistics):

        line = ",".join((chrm, pos, str(depth)))
        for pair in f_statistics.keys():
            if f_statistics[pair] == None:
                continue

            if f_statistics[pair][stat_id] == np.nan:
                value = 'NA'
            else:
                value = str(f_statistics[pair][stat_id])
            line += "," + value
        
        header = 'chrm,pos,total_depth,' + ','.join([left+"-"+right for left, right in f_statistics.keys()]) + "," + stat_id
        return (line, header)

    def phred2probablility(self, phred_score):
        if phred_score == "./.": 
            return phred_score
        else:
            phred_score = int(phred_score)
            #phred to probability = 10**(-phred/10)
            probability = 10**((phred_score*-1)/10.0)
        return probability

    def sample2population(self, sample_id):
        sample_pop = None
        for pop in self.populations.keys():
            if sample_id in self.populations[pop]:
                sample_pop = pop

        return sample_pop


    def slice_vcf(self, tabix_filename, chrm, start=None, stop=None):

        tabixfile = pysam.Tabixfile(tabix_filename)
        results = None
        # pysam throws exception when there are no SNPs in the range given.
        # this correction makes sure that some empty data is returned instead
        try:
            chunk = tabixfile.fetch(reference=chrm, start=start, end=stop)
        except:
            results = []

        if results != []:
            results = []
            for line in chunk:
                line = self.parse_vcf_line(line)
                results.append(line.copy())

        return results


    def calc_MAF(self, gt_likelihoods):
        pass


    def set_chrms(self, vcf_path):
        
        if vcf_path.endswith('gz'):
            vcf_file = gzip.open(vcf_path, 'rb')
        else:
            vcf_file = open(vcf_path,'rU')

        chrm_len_pairs = []
        for line in vcf_file:
            if line.startswith("##contig"):
                line_parts = line.strip("##contig=<ID=").strip(">\n").split(",")
                line_parts[-1] = int(line_parts[-1].strip("length="))
                chrm_len_pairs.append(line_parts)
            if line.startswith("#CHROM"):
                break
        
        vcf_file.close()
        self.chrm2length =  dict(chrm_len_pairs)

    def set_header(self, vcf_path):
        if vcf_path.endswith('gz'):
            vcf_file = gzip.open(vcf_path, 'rb')
        
        else:
            vcf_file = open(vcf_path,'rU')
        
        for line in vcf_file:
            if line.startswith("#CHROM"):
                self.header = line.strip("#").strip().split()
                self.__header_dict__ = OrderedDict([(item,None) for item in self.header])
                break

        vcf_file.close()

    def process_genotype(self, snp_call):
        """Given GT string convert to allele counts.

            e.g., '1/1' becomes (0,2)
                  '0/0' becomes (2,0)  
                  '0/1' becomes (1,1)                 

            To Do: Properly handle multiallelic GTs
        """

        if snp_call == "0/0":
            return (2,0)
        
        elif snp_call == "1/1":
            return (0,2)
        
        elif snp_call == '1/0' or \
               snp_call == '0/1':
            return (1,1)
        
        # skip multiallelic sites
        else:
            return (0,0)


    def get_major_allele(self, vcf_line):
        """Choose major allele from VCF line based on genotype counts."""

        samples = self.header[9:]
        GT_counts = (self.process_genotype(vcf_line[k]['GT']) for k in samples if vcf_line[k] != None)

        ref_count = sum(( gt[0] for gt in GT_counts))
        alt_count = sum(( gt[1] for gt in GT_counts))
        
        if alt_count > ref_count:
            return vcf_line['ALT']
        
        elif alt_count == ref_count:
            return choice((vcf_line['ALT'],vcf_line['REF']))
        
        else:
            return vcf_line['REF']

    def sequence_between_SNPS(self, fasta_handle, chrm, start, stop):
        """Returns sequence between SNPs. Requires an opened pysam faidxed fasta."""
        if start == None and stop == None:
            return None

        elif start == None and stop != None:
            start = 0
            stop = int(stop) - 1
            seq = fasta_handle.fetch(chrm,start, stop)
            return seq

        elif int(stop) - int(start) == 1:
            seq = ''
            return seq
  
        elif int(stop) - int(start) == 1:
            return ''
        
        else:
            start = int(start) + 1
            stop = int(stop) - 1
            seq = fasta_handle.fetch(chrm,start, stop)
            return seq

    def polarize_allele_counts(self, allele_counts, vcf_line):
        """For a set of allele counts use outgroup to determine REF and ALT alleles.

        The outgroup is required to be fixed at one allele. SNPS with variable outgroups are ignored.
        """

        # To Do: Check that outgroup is defined. Throw exception if is isn't. 

        # skip outgroup not fixed at one value
        if self.check_outgroup(allele_counts) == False:
            return None 

        # skip multi-allelic sites
        elif len(vcf_line['REF']) > 1 or len(vcf_line["ALT"]) > 1:
            return None                         

        else:
            # CALL BASE FOR OUTGROUP
            outgroup_allele = self.get_outgroup_base(allele_counts, vcf_line)
            
            # CALL MAJOR ALLELE (BASE) FOR INGROUP
            major_allele = self.get_ingroup_major_allele(allele_counts, vcf_line, outgroup_allele)

            # POLORIZE REF AND ALT FOR INGROUP
            if major_allele != vcf_line['REF']:
                ref, alt = ('ALT','REF')
            else:
                ref, alt = ('REF','ALT')

            polarized_calls = {"CHROM":vcf_line['CHROM'], 'POS':vcf_line['POS'], \
                               'MAJOR_ALLELE':major_allele, 'OUTGROUP_ALLELE':outgroup_allele, \
                               'REF':vcf_line['REF'], 'ALT':vcf_line["ALT"]}

            for count, pop in enumerate(self.populations.keys()):
                polarized_calls[pop] = {'REF':allele_counts[pop][ref], 'ALT':allele_counts[pop][alt]}

            return polarized_calls


    def count_alleles(self, vcf_line, polarize=False):
        """"Given a vcf file summarize allele counts per population."""

        allele_counts = {"CHROM":vcf_line['CHROM'], 'POS':vcf_line['POS']}

   
        for pop in self.populations.keys():
            
            alleles = {'REF':0, 'ALT':0}
            
            for sample in self.populations[pop]:
                if vcf_line[sample] != None:
                    ref, alt = self.process_genotype(vcf_line[sample]['GT'])
                    alleles['REF'] += ref
                    alleles['ALT'] += alt
            
            allele_counts[pop] = alleles.copy()

        if polarize == True:
            return self.polarize_allele_counts(allele_counts, vcf_line)

        else:
            return allele_counts
                

    def check_outgroup(self, row):
        outgroup = row["outgroups"]

        if outgroup['ALT'] == 0 and outgroup['REF'] > 0:
            return True

        elif outgroup['REF'] == 0 and outgroup['ALT'] > 0:
            return True

        else:
            return False


    def get_outgroup_base(self, row, raw_calls):
        pops = self.populations.keys()

        if row["outgroups"]['REF'] > row["outgroups"]['ALT']:
            outgroup_base = raw_calls['REF']
        
        elif sum(row['outgroups'].values()) == 0 : # e.g, == {'ALT': 0, 'REF': 0}
            ref = sum([row[pop]['REF'] for pop in pops])
            alt = sum([row[pop]['ALT'] for pop in pops])
            
            if ref > alt:
                outgroup_base = raw_calls['REF']
            else: 
                outgroup_base = raw_calls['ALT']
        
        else:
            outgroup_base = raw_calls['ALT']    
        
        return outgroup_base

    def get_ingroup_major_allele(self, row, raw_calls, outgroup_allele):
    
        pops = self.populations.keys()
        pops.remove('outgroups')

        ref_sum = sum([row[pop]['REF'] for pop in pops])
        alt_sum = sum([row[pop]['ALT'] for pop in pops])

        if ref_sum < alt_sum:
            return raw_calls['ALT']

        else:
            return raw_calls['REF']

    def make_triplet(self, base):
        return "-{0}-".format(base) 

    def count_alleles_in_vcf(self, vcf_slice):
        """"Fuction that returns rows of allele counts line by line of a given VCF."""

        # To Do: Add code to skip header. Currently assumes sliced vcf. 

        if vcf_slice == [] or vcf_slice == None:
            return None

        else:
            return [self.count_alleles(vcf_line, polarize=True) for vcf_line in vcf_slice]

    def sliding_window_dadi(self, args):

        # Using the data in the VCF header generate all slices
        print 'Generating Slices...'
        slices = generate_slices(args)
        
        # Open output file
        print 'Processing Slices...'
        fout = open(args.output,'w')

        # Write output header
        fout.write(','.join(create_header(pop_ids)) + "\n")

        for key_count, chrm in enumerate(slices.keys()):

            if args.region != [None]:
                chrm = args.region[0]
                #if key_count == 2: break
            
            for count, s in enumerate(slices[chrm]):

                # Break out of loop if loop proceeds beyond
                #   defined region (-L=1:1-5000 = 5000)
                if s[-1] > args.region[-1] and args.region[-1] != None: break

                region = [chrm] + list(s)

                # Setup Pairwise Dadi
                dadi_data = make_dadi_fs(args, region)
                if dadi_data == None: continue # skip empty calls

                dd, pop_ids = dadi_data
                projection_size = 10
                pairwise_fs  = dadi.Spectrum.from_data_dict(dd, pop_ids, [projection_size]*2)

                # Create final line, add Fst info
                final_line = region
                final_line += [pairwise_fs.Fst()]

                # Add in population level stats
                for pop in pop_ids:
                    fs = dadi.Spectrum.from_data_dict(dd, [pop], [projection_size])
                    final_line += [fs.Tajima_D(), fs.Watterson_theta(), fs.pi(), fs.S()]

                # write output
                final_line = [str(i) for i in final_line]
                fout.write(','.join(final_line) + "\n")

            # Don't process any more keys than necessary
            print 'Processed {0} of {1} slices from contig {2}'.format(count+1, len(slices[chrm]), chrm)
            if args.region != [None]: break

        fout.close()


    def generate_slices(self, args):

        if args.region != [None]:
            self.chrm2length = {args.region[0]:args.region[-1]}
        else:
            self.set_chrms(args.input)
        
        chrm_2_windows = self.chrm2length.fromkeys(self.chrm2length.keys(),None)

        for count, chrm in enumerate(self.chrm2length.keys()):

            length = self.chrm2length[chrm]
            window_size = args.window_size
            overlap = args.overlap

            # Skip contigs that are to short
            if length <= window_size: continue
            
            # Fit windows into remaining space
            if (length % window_size) > overlap:
                start = (length % window_size)/2
                stop = (length - window_size) - overlap/2

            # Prevent windows from invading remaining space 
            if (length % window_size) <= overlap:
                start = (length % window_size)/2
                stop = length - overlap*2
        
            if overlap == 0:
                starts = xrange(start, stop, window_size)
            else:
                starts = xrange(start, stop, overlap)
            
            stops = [i+window_size for i in starts]
            windows = zip(starts, stops)
            
            chrm_2_windows[chrm] = windows

        return chrm_2_windows


    def slice_2_allele_counts(self, tabix_filename, chrm, start=None, stop=None):
        vcf_slice = self.slice_vcf( tabix_filename, chrm, start, stop)
        processed_slice = self.count_alleles_in_vcf( vcf_slice)
        return processed_slice


    def vcf_slice_2_fstats(self, vcf_slice, projection_size = 10):

        pop_ids = self.populations.keys()
        pop_ids.remove('outgroups')

        processed_slice = self.count_alleles_in_vcf(vcf_slice)
        print self.calc_fstats(processed_slice)
        # calc_fstats
        #dd = self.make_dadi_fs(processed_slice)
        #fs = dadi.Spectrum.from_data_dict(dd, pop_ids, [projection_size]*len(pop_ids))
        #return fs.Fst()


    def make_dadi_fs( self, vcf_slice,):

        if vcf_slice == [] or vcf_slice == None:
            return None

        else:

            final_dadi = {}
            
            for row_count, row in enumerate(vcf_slice):
                
                if row == None: continue
                
                calls = {}
                for count, pop in enumerate(self.populations.keys()):
                    if pop == 'outgroups': continue
                    calls[pop] = (row[pop]['REF'], row[pop]['ALT'])

                row_id = "{0}_{1}".format(row['CHROM'],row['POS'])
            
                dadi_site = {'calls': calls,
                       'context': self.make_triplet(row['MAJOR_ALLELE']),
                       'outgroup_context': self.make_triplet(row['OUTGROUP_ALLELE']),
                       'outgroup_allele': row['OUTGROUP_ALLELE'],
                       'segregating': (row['REF'], row['ALT'])
                       }

                final_dadi[row_id] = dadi_site

            return final_dadi



    def SNPwise_fstats(self, vcf, vcf_slice):

        #vcf, vcf_slice = vcf_slice
        results = []
        for vcf_line in vcf_slice:
            allele_counts = vcf.count_alleles(vcf_line, polarize=False)
            try:
                results.append(vcf.calc_fstats(allele_counts))
            except:
                continue

        return results


