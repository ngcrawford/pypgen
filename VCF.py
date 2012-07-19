import gzip
import pysam
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

    def parse_individual_snps(self, vcf, fout, stat_id):

        vcf = open(vcf,'rU')
        
        if fout == None:
            fout = sys.stdout
        else:
            fout = open(fout,'w')

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
  
                allele_counts = self.calc_allele_counts(vcf_line_dict)
                f_statistics = self.calc_fstats(allele_counts)
                info = self.parse_info(vcf_line_dict["INFO"])
                formated_data, header = self.__write_to_outfiles__(vcf_line_dict["CHROM"],\
                                    vcf_line_dict["POS"],info["DP"], stat_id, f_statistics)

                if line_count == 1:
                    fout.write(header+"\n")
                    fout.write(formated_data+"\n")
                else:
                    fout.write(formated_data+"\n")

            if line_count >= 0:
                line_count += 1

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
            # self.sample_format_dict = dict([(item,None) for item in self.sample_format])      
     
        for count, item in enumerate(vcf_line_dict):
            
            if count >= 9:
                genotype = vcf_line_dict[item]
                
                if genotype == "./.":
                    vcf_line_dict[item] = None

                else:
                    genotype = dict(zip(self.sample_format,genotype.split(":")))
                    vcf_line_dict[item] = genotype

        return vcf_line_dict

    def calc_allele_counts(self, vcf_line_dict):

        allele_counts = self.populations.fromkeys(self.populations.keys(),None)

        for population in self.populations.keys():

            allele_format_dict = {0:0,1:0,2:0,3:0,4:0}   # create dict to prevent pointer issues
            allele_counts[population] = allele_format_dict

            for sample_id in self.populations[population]:
                

                if vcf_line_dict[sample_id] != None:
                    
                    genotype = vcf_line_dict[sample_id]
                    genotype = genotype["GT"].split("/")

                    if genotype == [".","."]: continue

                    genotype = [int(item) for item in genotype]
                    
                    for allele in genotype:
                        allele_counts[population][allele] += 1 
        
        return allele_counts

   
    def calc_fstats(self, allele_counts):

        #test_pair = {'So_Riviere_Goyaves': np.array([ 0.0, 1.0, 0.0, 0.0]), 'Plage_de_Viard': np.array([ 1.0, 0.0, 0.0, 0.0]),}
        
        # CALCULATE ALLELE FREQUENCIES
        allele_freqs_dict = self.populations.fromkeys(self.populations.keys(),None)
        for population in allele_counts.keys():
            counts =  allele_counts[population].values()
            freqs =  counts/np.sum(counts,dtype=float)
            allele_freqs_dict[population] = freqs

        allele_freqs = allele_freqs_dict.values()

        # CACULATE PAIRWISE F-STATISTICS
        pairwise_results = {}
        for population_pair in combinations(self.populations.keys(),2):
            
            pop1, pop2 =  population_pair
            Ns = [sum(allele_counts[pop].values()) for pop in [pop1, pop2]]

            if 0 in Ns: 
                values = [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]
                values_dict = dict(zip(['Hs_est', 'Ht_est', 'Gst_est', 'G_prime_st_est', 'G_double_prime_st_est', 'D_est'],values))
                pairwise_results[population_pair] = values_dict
                continue
            
            else:

                pop1 = allele_freqs_dict[pop1]
                pop2 = allele_freqs_dict[pop2]
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


    def slice_vcf(self, tabix_filename, pops, chrm, start=None, stop=None):

        tabixfile = pysam.Tabixfile(tabix_filename)
        try:
            chunk = tabixfile.fetch(reference=chrm, start=start, end=stop)
        except:
            return None 

        else:           
            pop_ids = self.populations.keys()
            lines = []

            for line in chunk:
                line =  self.parse_vcf_line(line)
                lines.append(line.copy())

            return lines


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

            if line.startswith("#CHROM"): break

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

