#!/usr/bin/env python
# encoding: utf-8
"""
pypgen.py

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

The author may be contacted at ngcrawford@gmail.comm
"""

import os
import re
import csv
import sys
import copy
import unittest
import argparse
import pandas
from la import *
from pylab import *
import multiprocessing
from functools import partial



def get_args():
    """Parse sys.argv"""
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input-file', required=True, 
        help='The input file. Either genepop or arlequin formated.')
    args = parser.parse_args()
    return args


class population(object):
    """docstring for population"""
    def __init__(self, genotypes, name=None, loci=None, individuals=None,\
            missing_data=('000000', '0000', 'NA', '00', '0', '??', '?', 'BADDNA')):
        super(population, self).__init__()
       
        def label_genotypes(self):
            """Convert an unlabeled array of genotypes to a larry"""
            
            if self.individuals == None:    # give numeric names (0-some #) if individuals are unnamed.
                individuals = range(0,len(self.genotypes))
                labeled_genotypes = larry(self.genotypes, [individuals, self.loci])
            else:
                labeled_genotypes = larry(self.genotypes, [self.individuals, self.loci])
            return labeled_genotypes
        
        # population data
        self.missing_data = missing_data
        self.name = name
        self.loci = loci
        self.individuals = individuals
        self.genotypes = array(genotypes)
        self.genotypes = label_genotypes(self)   # store genotypes as labeled array (larry)
    
    
    def __empty_matrix__(self):
        """ creat empty matrix for counts or frequencies"""
        unique_alleles = self.unique_alleles()
        ua_len, l_len = (len(unique_alleles), len(self.loci))
        empty_data = zeros((ua_len ,l_len))
        allele_counts = larry(empty_data,[unique_alleles,self.loci])
        return allele_counts
    
    def exp_het(self):
        "Calculate Expected Heterozygosity"
        Hexp_dict = {}
        for locus in self.loci:
            allele_freqs = self.allele_freqs().lix[:,[locus]]
            allele_freqs_squared = pow(allele_freqs,2)              # square freqs
            Hexp_dict[locus] = 1-(allele_freqs_squared.sum())
        return Hexp_dict.copy()
    
    def individuals_per_locus(self):
        """Count the number of individuals per locus."""
        genotypes = self.genotypes
        loci = genotypes.getlabel(1)
        empty_array = zeros((len(loci)))
        storage_larry = larry(empty_array,[loci])
        
        for locus in loci:
            geno_at_locus = genotypes.lix[:,[locus]]
            individual_count = 0
            for genotype in geno_at_locus:
                if genotype in self.missing_data: continue # ignore missing data
                else: individual_count += 1
            
            storage_larry.set([locus], individual_count)
        return storage_larry
        
    def unique_alleles(self):
        """Get list of unique alleles in population of genotypes"""
        alleles = []
        for individual in self.genotypes:
            for genotype in individual:
                if genotype in self.missing_data: continue # ignore missing data
                    
                left = genotype[:len(genotype)/2]   # split allele pairs into left and right pieces
                right = genotype[len(genotype)/2:]
                alleles.append(left)
                alleles.append(right)

        unique_alleles = list(set(alleles))
        return unique_alleles
        
    def allele_counts(self):
        """Convert set of raw genotypes (larry) to a matrix of allele counts"""
        # set up empty labeled set of data
        allele_counts = self.__empty_matrix__()
            
        # added counts to alleles
        for locus in self.loci:
            for genotype in self.genotypes.lix[:,[locus]]:
                if genotype in self.missing_data: continue # ignore missing data
                
                left = genotype[:len(genotype)/2]   # split allele pairs into left and right pieces
                right = genotype[len(genotype)/2:]
                
                current_count = allele_counts.get([left,locus]) # update allele_counts
                current_count += 1
                allele_counts.set([left,locus],current_count)
                
                current_count = allele_counts.get([right,locus]) # update allele_counts
                current_count += 1
                allele_counts.set([right,locus],current_count)
                
        return allele_counts
    
    def allele_freqs(self):
        """Convert matrix of allele counts to allele frequencies"""
         # set up empty labeled set of data
        allele_freqs = self.__empty_matrix__()
        
        for locus in self.loci:
            locus_counts = self.allele_counts().lix[:,[locus]]
            total_count = locus_counts.sum()
            
            for count, allele_count in enumerate(locus_counts):
                freq = allele_count/total_count
                locus_labels = locus_counts.getlabel(axis=0)        # this is a bit hacked, could be improved
                allele_freqs.set([locus_labels[count],locus],freq)
        return allele_freqs
        
    def n(self):
        """Returns a dictionary of number of individuals with 
        at a particular locus with genotypes. Individuals with missing
        data are ignored.
        
        Examples
        --------
        
        more here...
        
        {'Locus 1': 10, 'Locus 2': 10}
        """
        n = {}
        for locus in self.loci:
            locus_len = 0
            for genotype in self.genotypes.lix[:,[locus]]:
                if genotype in self.missing_data: continue # ignore missing data
                else:
                    locus_len += 1
            n[locus] = locus_len
        n = n.copy()
        return n

class populations(list):
    """Class to operate """
    def __init__(self,):
        super(populations, self).__init__()
        self.pops = []
    
    def __iter__(self):
        """Make the pops in populations an interable"""
        for pop in self.pops:
            yield pop
    
    def __getitem__(self, key):
        return self.pops[key]
            
    def __allele_freqs_3D_larry__(self):
        "returns populations as labeled 'larry.'"
        pop_names = self.pop_names()
        pop_loci = self.pops[0].loci
        unique_alleles = self.unique_alleles()
        empty_data = zeros((len(unique_alleles), len(pop_loci)))
        complete_alleles = larry(empty_data.copy(),[unique_alleles, pop_loci])
        data3d = []
        for count, pop in enumerate(self.pops):
            pop.allele_freqs()
            alleles_data = complete_alleles.merge(pop.allele_freqs(), update=True)
            data3d.append(alleles_data.copyx())
        data3d = array(data3d)
        data3d = larry(data3d.copy(), [pop_names,unique_alleles,pop_loci])
        return data3d
    
    def __empty_pop_by_loci_larry__(self):
        """Make Empty Array for Internal Use.
        should probably be __empty_pop_by_loci_larry__"""
        pop_names = self.pop_names()
        empty_data = zeros((len(pop_names), len(self.pops[0].loci)))
        labeled_empty_larry = larry(empty_data.copy(), [pop_names, self.pops[0].loci])
        return labeled_empty_larry
        
    def pop_names(self):
        pop_names = []
        for pop in self.pops:
            pop_names.append(pop.name)
        return pop_names
    
    def unique_alleles(self):
        "Returns all allele names in all populations removing duplicates"
        allele_names = []
        for pop in self.pops:
            for allele in pop.unique_alleles():
                allele_names.append(allele)
        
        allele_names = list(set(allele_names))          # remove duplicate names
        return allele_names
        
    def allele_counts(self):
        """Count the number of loci in each population 
        accounting for missing data."""

        n_alleles = self.__empty_pop_by_loci_larry__()
        for pop in self.pops:
            n_dict = pop.n()
            for locus in self.pops[0].loci:
                n = n_dict[locus]
                n_alleles.set((pop.name, locus), n)     
        return n_alleles
    
    def proportion_variable_loci(self):
        """Unfinished: prints out the population name and the propotion of variable loci"""
        
        for pop in self.pops:
            allele_counts = pop.allele_counts()
            variable_loci = 0.0
            monomorphoic_loci = 0.0
            alleles_in_pop = 0
            for pos in range(allele_counts.shape[-1]):
                if len(set(allele_counts[:,pos])) == 2:
                    monomorphoic_loci += 1
                else:
                    variable_loci += 1
                alleles_in_pop +=  len(set(allele_counts[:,pos]))
            print pop.name, variable_loci/len(pop.loci), alleles_in_pop
        
    def loci_harmonic_means(self):
        """calculates harmonic mean from list of integers"""
        
        loci = self.pops[0].loci
        empty_array = zeros((len(loci)))
        fractional_allele_counts = larry(empty_array,[loci])
        
        n = 0
        for pop in self.pops:
            inds_per_locs = pop.individuals_per_locus()
            fractional_counts = 1/inds_per_locs
            fractional_allele_counts += fractional_counts
            n += 1
        loci_harmonic_means = n/fractional_allele_counts

        return loci_harmonic_means

    def harmonic_mean_chao(self, values):
        """Calculates the harmonic mean following the method suggested by
        Anne Chao. The formula is: 1/[(1/A)+var(D)(1/A)**3]. Used for 
        calculating multilocus Dest."""
        A = values.mean()
        varD = values.var()
        harmonic_mean_chao = 1/((1/A)+(varD)*pow((1/A),3))
        return harmonic_mean_chao
    
    def n(self):
        n = float(len(self.pop_names()))
        return n
    
    
    def Hs(self, by_pop=False):
        """Calculate Hs (mean within-subpopulation heterozygosity, Nei and Chesser 1983)"""
        # setup names, empty array, and final storage matrix
        heterozygosity = self.__empty_pop_by_loci_larry__()
        
        # loop through population calculating heterozygosity for each locus.
        # Each het. then gets 'set', really appended, to the final
        # storage array
        for pop in self.pops:
            pop_heterozygosities = pop.exp_het()
            for locus in pop_heterozygosities.keys():
                current_het = pop_heterozygosities[locus]
                heterozygosity.set((pop.name, locus), current_het)
        
        # return the mean heterozygosity for each locus.
        if by_pop == True:
            print heterozygosity.std(axis=1)
            return heterozygosity.mean(axis=1)
        else:    
            return heterozygosity.mean(axis=0)
               
    def Ht(self):
        """Calculate Ht the heterozygosity of the pooled subpopulations (Nei and Chesser 1983)"""
        # for the same allele in each locus sum the frequencies from each population and square
        # then sum all those frequencies and subtract from 1.0
        
        empty_data =  zeros((len(self.pops[0].loci)))
        Ht_larry = larry(empty_data.copy(), [self.pops[0].loci])

        for locus in self.pops[0].loci:
            allele_names = self.unique_alleles()
            pop_names = self.pop_names()
            empty_data = zeros((len(allele_names), len(pop_names)))
            labeled_empty_larry = larry(empty_data.copy(), [allele_names, pop_names])
            
            for pop in self.pops:
                pop_allele_freqs = pop.allele_freqs()
                locus_freqs = pop_allele_freqs.lix[:,[locus]]
                allele_names = locus_freqs.getlabel(axis=0)
            
            Ht = 1.0-(pow(labeled_empty_larry.mean(axis=1),2).sum())
            Ht_larry.set([locus],Ht)
            
        return Ht_larry
       
    def Hst(self):
        Ht = self.Ht()
        Hs = self.Hs()
        Hst = (Ht-Hs)/(1.0-Hs)
        return Hst
    
    def delta_s(self):
        Hs = self.Hs()
        delta_s = pow((1.0-Hs),-1.0)
        return delta_s
    
    def delta_t(self):
        Ht = self.Ht()
        delta_t = pow((1.0-Ht),-1.0)
        return delta_t
    
    def Dst(self):
        """Nei's Dst (absolute differentiation)"""
        Ht = self.Ht()
        Hs = self.Hs()
        Dst = Ht-Hs
        return Dst
    
    def Gst(self):
        """Nei's Gst (relative differentiation)"""
        Ht = self.Ht()
        Dst = self.Dst()
        Gst = (Dst)/Ht     
        return Gst
    
    def D(self):
        """Jost's D uncorrected for sample size."""
        Ht = self.Ht()
        Hs = self.Hs()
        n = self.n()
        D = ((Ht-Hs)/(1.0-Hs))*(n/(n-1.0))
        return D
    
    # DIVERSITY ESTIMATORS
    def Hs_prime_est(self):
        """Calculate corrected Hs: the mean within-subpopulation heterozygosity (Nei and Chesser 1983)."""
        n = self.n()
        allele_freqs = self.__allele_freqs_3D_larry__()
        Hj = 1-(allele_freqs.power(2).sum(axis=1))
        Hs_prime_est = (1/n)*(Hj.sum(axis=0))
        return Hs_prime_est
    
    def Ht_prime_est(self):
        """Calculate corrected Ht: the heterozygosity of the pooled subpopulations (Nei and Chesser 1983)"""
        n = self.n()
        allele_freqs = self.__allele_freqs_3D_larry__()
        inner = ((1/n)*allele_freqs.sum(axis=0)).power(2)
        Ht_prime_est = 1-inner.sum(axis=0)
        return Ht_prime_est
    
    def Hs_est(self):
        """ Basic Equation: ((2*N_harmonic)/(2*N_harmonic-1))*Hs"""
        locus_harmonic_means = self.loci_harmonic_means()
        Hs = self.Hs_prime_est()
        Hs_est = ((2.0*locus_harmonic_means)/(2.0*locus_harmonic_means-1.0))*Hs
        return Hs_est
    
    def Ht_est(self):
        """Basic Equation: Ht+Hs_est/(2*N_harmonic*n)"""
        locus_harmonic_means = self.loci_harmonic_means()
        Hs_est = self.Hs_est()
        Ht = self.Ht_prime_est()
        n = len(self.pop_names())
        Ht_est = Ht+Hs_est/(2.0*locus_harmonic_means*n)
        return Ht_est
    
    def Gst_est(self):
        Ht_est = self.Ht_est()
        Hs_est = self.Hs_est()
        Gst_est = (Ht_est-Hs_est)/Ht_est
        return Gst_est
    
    def G_prime_st_est(self):
        Ht_est = self.Ht_est()
        Hs_est = self.Hs_est()
        G_est = self.Gst_est()
        n = self.n()
        G_prime_st = (G_est*(n-1.0+Hs_est))/((n-1.0)*(1.0-Hs_est))
        return G_prime_st
    
    def G_double_prime_st_est(self):
        """G''st = k*(HT-HS)/((k*HT-HS)*(1-HS)"""
        Ht_est = self.Ht_est()
        Hs_est = self.Hs_est()
        n = self.n()
        G_double_prime_st_est = n*(Ht_est-Hs_est)/((n*Ht_est-Hs_est)*(1-Hs_est))
        return  G_double_prime_st_est
    
    def D_est(self):
        Ht_est = self.Ht_est()
        Hs_est = self.Hs_est()
        n = self.n()
        D_est = ((Ht_est-Hs_est)/(1.0-Hs_est))*(n/(n-1))
        return D_est
    
    def D_est_chao(self):


        def a(self):
            for locus in self.pops[0].loci:
                for pop in self.pops:
                    pop.lix[:,[locus]]
                    
                
            pass
        
        def b(self):
            pass
        
        
        empty_data =  zeros((len(self.pops[0].loci)))
        storage_loci_larry = larry(empty_data.copy(), [self.pops[0].loci])
        pass
        
    
    # MULTILOCUS FUNCTIONS
    def multilocusGst_est(self):
        """Averages across loci before calculating Gst."""
        Ht_est = self.Ht_est().mean()
        Hs_est = self.Hs_est().mean()
        multilocusGst_est = (Ht_est-Hs_est)/Ht_est
        return multilocusGst_est
    
    def multilocusG_prime_st_est(self):
        Ht_est = self.Ht_est().mean()
        Hs_est = self.Hs_est().mean()
        G_est = self.multilocusGst_est()
        n = self.n()
        G_prime_st = (G_est*(n-1.0+Hs_est))/((n-1.0)*(1.0-Hs_est))
        return G_prime_st
    
    def multilocusG_double_prime_st_est(self):
        Ht_est = self.Ht_est().mean() # take mean across loci
        Hs_est = self.Hs_est().mean()
        G_est = self.multilocusGst_est()
        n = self.n()
        G_double_prime_st_est = n*(Ht_est-Hs_est)/((n*Ht_est-Hs_est)*(1.0-Hs_est))
        return G_double_prime_st_est
    
    def multilocusD_est(self):
        D_est = self.D_est()
        multilocusD_est = self.harmonic_mean_chao(D_est)
        return multilocusD_est
    
    def all_estimators(self):
        """Calculate all the estimators simultaneously (e.g., Gst, G'st, G''st and Dest)
        This is substantially faster than calculating each one independantly."""
        n = self.n()
        Ht_est = self.Ht_est()
        Hs_est = self.Hs_est()
        
        # calculate estimators. Saving values as narrays using .copyx()
        Gst_est = ((Ht_est-Hs_est)/Ht_est).copyx()
        G_prime_st_est = ((Gst_est*(n-1.0+Hs_est))/((n-1.0)*(1.0-Hs_est))).copyx()
        G_double_prime_st_est = (n*(Ht_est-Hs_est)/((n*Ht_est-Hs_est)*(1-Hs_est))).copyx()
        D_est = (((Ht_est-Hs_est)/(1.0-Hs_est))*(n/(n-1))).copyx()
        
        # create final larry with appropriate labels
        estimator_array = larry(array([Gst_est, G_prime_st_est, G_double_prime_st_est, D_est]),
                                [["Gst-est", "G'st-est", "G''st-est","D-est"], Ht_est.getlabel(axis=0)])
        return estimator_array

    
    def all_multilocus_estimators(self, n=None, Ht_est=None, Hs_est=None):
        """Calculate all the multilocus estimators simultaneously (e.g., Gst, G'st, G''st and Dest)
        This is substantially faster than calculating each one independantly."""
        
        if n == None and Ht_est == None and Hs_est == None:
            n = self.n()
            Ht_est = self.Ht_est()
            Hs_est = self.Hs_est()
            Ht_est_mean = Ht_est.mean()
            Hs_est_mean = Hs_est.mean()
        
        # Calculate multilocus estimators using the means of the Ht and Hs estimators 
        #   Note that Dest is estimated using the harmonic mean chao function
        multilocus_Gst_est = (Ht_est_mean-Hs_est_mean)/Ht_est_mean
        multilocus_G_prime_st = (multilocus_Gst_est*(n-1.0+Hs_est_mean))/\
                                    ((n-1.0)*(1.0-Hs_est_mean))
        multilocus_G_double_prime_st_est = n*(Ht_est_mean-Hs_est_mean)/((n*Ht_est_mean-Hs_est_mean)*(1.0-Hs_est_mean))
        multilocus_D_est = self.harmonic_mean_chao(((Ht_est-Hs_est)/(1.0-Hs_est))*(n/(n-1)))
        
        multilocus_estimators_array = larry(array([multilocus_Gst_est,multilocus_G_prime_st,multilocus_G_double_prime_st_est,multilocus_D_est]),
                                            [["multilocus Gst-est","multilocus G'st-est","multilocus G''st-est","multilocus D-est"]])

        return multilocus_estimators_array
    
    def __randomize_pop_pairs__(self, z, y):
 
        z_g = z.genotypes.copyx()
        y_g = y.genotypes.copyx()
  
        merged_data = np.row_stack((z_g,y_g))

        shuffled_loci = array([np.random.permutation(item) \
                        for item in merged_data.transpose()]).transpose() # shuffle columns

        z_shuffled = shuffled_loci[:z_g.shape[0]]
        y_shuffled = shuffled_loci[z_g.shape[0]:]

        z_shuffled_pop = population(z_shuffled , loci=z.loci, name=z.name)
        y_shuffled_pop = population(y_shuffled , loci=y.loci, name=y.name)


        return (z_shuffled_pop,y_shuffled_pop)

    def calculate_pairwise_estimators(self, estimator='D_est', randomize_samples=False):
        """rename"""
        pop_names = len(self.pop_names())
        storage_array = larry(zeros((pop_names,pop_names)), [self.pop_names(), self.pop_names()])
        
        for count, pop_list in enumerate(self.pops):
            count = count + 1 # correct offset
            for left, right in zip(self.pops[count:], self.pops[:-count]):
                
                # setup population class and add paired populations
                paired_pops = populations()

                if randomize_samples == True:
                    right, left = self.__randomize_pop_pairs__(right, left)

                paired_pops.append(right)
                paired_pops.append(left)
                
                # apply appropriate multilocus estimator
                if estimator == 'Gst_est':
                    value = paired_pops.multilocusGst_est()
                if estimator == 'G_prime_st_est':
                    value = paired_pops.multilocusG_prime_st_est()
                if estimator == 'G_double_prime_st_est':
                    value = paired_pops.multilocusG_double_prime_st_est()
                if estimator == 'D_est':
                    value = paired_pops.multilocusD_est()
                
                # update storage array
                storage_array.set(paired_pops.pop_names(),value)

        return storage_array
     
    def pairwise_estimators(self, estimator='D_est'):
        print estimator
        # get data, header, and side label
        data = self.calculate_pairwise_estimators(estimator)
        header = data.getlabel(0)
        side = data.getlabel(1)
        header = [estimator] + list(header) # add estimator labeled cell to header 
        
        # set up csv writer and write header
        filename = '/Users/ngcrawford/%s.csv' % (estimator)
        outfile = open(filename, 'w')
        data_writer = csv.writer(outfile)
        data_writer.writerow(header)
        
        # write data
        for count, row in enumerate(data):
            side_item = side[count]
            row = [side_item] + list(row)
            data_writer.writerow(row)
    
   
    def make_bootstraps(self, replicates=10, start=0):
        
        bootreps = []
        for rep in range(0,replicates):
            booted_pops = populations()

            for count, pop in enumerate(self.pops):
                sample_size = pop.genotypes.shape[0]
                choices = list(np.random.random_integers(0, sample_size-1, sample_size))
                
                # bootstrap individuals within each population
                booted_genotypes = pop.genotypes.copyx()[choices,:]
                new_pop = population(booted_genotypes, loci=pop.loci, name=pop.name) 
                booted_pops.append(new_pop)
            
            bootreps.append(booted_pops)
        
        return bootreps


    def bootstrap_all_multilocus_estimators(self):
        bootstraps = self.make_bootstraps()
        for pop in bootstraps:
            print pop.all_multilocus_estimators()

    def bootstrap_pairwise_estimators(self):
        bootstraps = self.make_bootstraps(replicates=200)
        pairwise_estimators = []
        
        for count, pop in enumerate(bootstraps):
            print count
            pairwise_estimators.append(pop.calculate_pairwise_estimators("Gst_est"))

        return pairwise_estimators

    def bootstrap_all_estimators(self):
        bootstraps = self.make_bootstraps(replicates=200)
        all_estimators = []
        for count, pops in enumerate(bootstraps):
            print count
            # print pops.all_estimators()
            all_estimators.append(pops.all_estimators())
        
        return all_estimators
        
    def write_genepop():
        """Write demes class to GenePop file"""
        pass
    
    def write_Rstat():
        """Write output for input into DEMEtics"""
    
    def read_arlequin():
        """Read in arlequin file into demes class"""
        pass
    
    def write_arlequin():
        """Write demes class to Arlequin  file"""
        pass
    
    def write_CSV():
        """Write demes class to comma/tab/other delimited file"""
        pass
    
    def write_genalex():
        """Write demes class to CSV file suitable for input into genalex"""
        pass
    
    def append(self, pop):
        # test if population
        self.pops.append(pop)
        # added code to update matrix too...
    def from_list(self,pops):
        self.pops = pops


# THERES A BUG IN PYTHON THAT PREVENTS A 'POOL'
# FROM CORRECTLY IMPORTING A CLASS FUNCTION
# PUTTING THE FUNCTION BEING 'MAPPED' CORRECTS THIS,
# BUT IT IS CLUDGY.

def wrap_calculate_pairwise_estimators(args):
    """Wrapper function that takes a tuple of argments
       and is therefore appopriate for mapping in a 
       multiprocessing pool."""
       
    demes, count, estimator = args
    print 'running replicate {0}'.format(count) 
    d = demes.calculate_pairwise_estimators(estimator,\
        randomize_samples=True)
    return d

class MonteCarlo(object):
    """Class for simulating null distributions and estimating p-values
        from pairs of populations.
        
        Requires Python 2.7 and for "pandas" to be installed
        
        """
    
    def __init__(self,):
        super(MonteCarlo, self).__init__()

    def simulate(self, demes, replicates=4, estimator="Gst_est"):
        """Simulate pairwise multi-locus estimators."""
        
        p = multiprocessing.Pool(multiprocessing.cpu_count())

        result = p.map(wrap_calculate_pairwise_estimators,\
            [(demes, count, estimator) for count in range(0,replicates)])

        return result

    def p_values(self, data, demes, estimator="D_est"):
        
        # get values from larrys
        ready_4_panda = []
        for count, datum in enumerate(data):
            ready_4_panda.append(datum.copyx())
        ready_4_panda = array(ready_4_panda)

        # get population names
        pop_names = data[0].label[0]

        # create panda (panal data array)
        final = pandas.Panel(ready_4_panda,\
            major_axis= pop_names,\
            minor_axis= pop_names)
            
        # cacluate observed values and convert to panda
        obs = demes.calculate_pairwise_estimators(estimator=estimator)
        obs = pandas.DataFrame(obs.copyx(),
            index=pop_names, 
            columns=pop_names)

        # Do pairwise p_value calculation
        for count, pop_list in enumerate(pop_names):
            count = count + 1 # correct offset
            for left, right in zip(pop_names[count:], pop_names[:-count]):
                null_dist = final.major_xs(right).xs(left)
                real_value = obs[left][right]
                p_value = 1.0 - np.searchsorted(null_dist,real_value)\
                                  /float(null_dist.shape[0])
                obs[right][left]= p_value
                print p_value
        
        print obs
        # store = pandas.HDFStore('store.h5')
        # store['final'] = final
        
def parse_genepop(lines):
    """process genpop lines into multidimentional array"""
    
    header = lines[0].strip()           # header info
    loci_names = []
    all_loci = []
    loci_list = []
    population_names = []
    pops_flag = 0                       # zero indicates in population names
    population_name = ''                # use to find first line of a population
    loci_on_multiple_lines = True
    line_counter = 0                    # tracks index of current line
    
    # REGEX stuff
    pattern = re.compile('pop', re.IGNORECASE)
    loci_punct = re.compile(',(?:\s*)|\s*')

    # GenePop Format Test
    if pattern.match(lines[2]):                     # assumes multiple populations...
        loci_names = loci_punct.split(lines[1])
        loci_on_multiple_lines = False

    for line in lines[1:]:                          # skip header line
        line = line.strip()                         # clean whitespace
        if len(line) == 0: continue                 # skip blank lines
        line_match = pattern.match(line)

        # DO STUFF WITH POPULATIONS
        if pops_flag >= 1:
            if line_match == None:
                line_parts = line.split(',')                        # split line with populations
                loci = line_parts[1].strip()                        # id loci
                loci_list.append(loci.split())                      # add loci to loci list

        # FLAG POPULATION AND UPDATE STORAGE LIST (All_Loci)
        if line_match:              # Flag Pops
            population_name = lines[line_counter+2].split()[0]      # get current pop name
            population_name = population_name.strip(',')
            population_names.append(population_name)
            pops_flag += 1
            if pops_flag >= 2:
                all_loci.append(loci_list)
            loci_list = []

        # PUT POPULATION NAMES IN A LIST
        if pops_flag == 0:              # get population names
            if loci_on_multiple_lines == True:
                loci_names.append(line)     # make list of population names
        line_counter += 1
    all_loci.append(loci_list)          # add last set of loci to all_loci

    # create populations class
    demes = populations()
    loci_names, population_names, all_loci
    for count, pop in enumerate(all_loci):
        new_pop = population(pop, loci=loci_names, name=population_names[count])  
        demes.append(new_pop)
    return demes




# Unit Tests
class InputFileTest(unittest.TestCase):
    genpopformat = """Title line: delete this example..
    Locus 1
    Locus 2
    Pop
    A , 001001 002002
    A , 001001 002002
    A , 002002 002002
    A , 002002 002002
    A , 002002 002002
    A , 002002 001001
    A , 002002 001001
    A , 002002 001001
    A , 002002 001001
    A , 002002 001001
    Pop
    B , 002002 002002
    B , 002002 002002
    B , 001001 002002
    B , 001001 002002
    B , 001001 002002
    B , 001001 001001
    B , 001001 001001
    B , 001001 001001
    B , 001001 001001
    B , 001001 001001
    """
    
    def testReadGenPop(self):
        genpopformat =  self.genpopformat
        lines = []
        for line in genpopformat.split('\n'):
            line = line.strip()
            lines.append(line)
        data = parse_genepop(lines)
        
        
        
class PopulationTests(unittest.TestCase):
    # test genotypes
    popA = [['001001', '002002'],
            ['001001', '002002'],
            ['002002', '002002'],
            ['002002', '002002'],
            ['002002', '002002'],
            ['002002', '001001'],
            ['002002', '001001'],
            ['002002', '001001'],
            ['002002', '001001'],
            ['002002', '001001']]

    def testUnique_alleles(self):
        pop = population(self.popA, loci=['Locus 1', 'Locus 2'])
        unique_alleles = pop.allele_counts()
        testvalues = larry.fromtuples([('002', 'Locus 1', 16.0), 
                                        ('002', 'Locus 2', 10.0), 
                                        ('001', 'Locus 1', 4.0), 
                                        ('001', 'Locus 2', 10.0)])
        self.assertEqual(unique_alleles,testvalues)
    
    def testAlleleCounts(self):
        pop = population(self.popA, loci=['Locus 1', 'Locus 2'])
        allelecounts = pop.allele_counts()
        testvalues = larry.fromtuples([('002', 'Locus 1', 8.0), 
                                        ('002', 'Locus 2', 5.0),
                                        ('001', 'Locus 1', 2.0),
                                        ('001', 'Locus 2', 5.0)])
        self.assertEqual(allelecounts,testvalues) # not working
        
    def testAlleleFreqs(self):
        pop = population(self.popA, loci=['Locus 1', 'Locus 2'])
        allelefreqs = pop.allele_freqs()
        testvalues = larry.fromtuples([('002', 'Locus 1', 0.80000000000000004),
                                        ('002', 'Locus 2', 0.5), 
                                        ('001', 'Locus 1', 0.20000000000000001),
                                        ('001', 'Locus 2', 0.5)])
        result = allelefreqs == testvalues
        for item in result.tolist()[0]:
            self.assertTrue(item,'Allele frequencies do not match')
    
    def testExpectedHeterozygosity(self):
        pop = population(self.popA, loci=['Locus 1', 'Locus 2'])
        Hexp = pop.exp_het()
        testvalues = {'Locus 1': 0.31999999999999984, 'Locus 2': 0.5}
        self.assertEqual(Hexp,testvalues) # not working
    
    def testPopSizeCalc(self):
        pop = population(self.popA, loci=['Locus 1', 'Locus 2'])
        testvalues = {'Locus 1': 10, 'Locus 2': 10}
        self.assertEqual(pop.n(), {'Locus 1': 10, 'Locus 2': 10})

class PopulationsTests(unittest.TestCase):
    
    #        Locus 1    Locus 2
    popA = [['001001', '002002'],
            ['001001', '002002'],
            ['002002', '002002'],
            ['002002', '002002'],
            ['002002', '002002'],
            ['002002', '001001'],
            ['002002', '001001'],
            ['002002', '001001'],
            ['002002', '001001'],
            ['002002', '001001']]

    popB = [['002002', '002002'],
            ['002002', '002002'],
            ['001001', '002002'],
            ['001001', '002002'],
            ['001001', '002002'],
            ['001001', '001001'],
            ['001001', '001001'],
            ['001001', '001001'],
            ['001001', '001001'],
            ['001001', '001001']]
            
    def make_test_demes(self):
        """setups up basic test deme"""
        pop1 = population(self.popA, loci=['Locus 1', 'Locus 2'], name = 'PopA')
        pop2 = population(self.popB, loci=['Locus 1', 'Locus 2'], name = 'PopB')
        testdemes = populations()
        testdemes.append(pop1)
        testdemes.append(pop2)
        return testdemes

    def testHs(self):
        testdemes = self.make_test_demes()
        Hs = testdemes.Hs()
        testvalues = larry.fromtuples([('Locus 1', 0.31999999999999984), 
                                       ('Locus 2', 0.5)])
        self.assertEqual(Hs, testvalues, 'Incorrect Hs values')
    
    def testHs_prime_est(self):
        testdemes = self.make_test_demes()
        Hs_prime_est = testdemes.Hs_prime_est()
        testvalues = larry.fromtuples([('Locus 1', 0.31999999999999984), 
                                      ('Locus 2', 0.5)])
        self.assertEqual(Hs_prime_est, testvalues, "Incorrect Hs'-est values")                       
        
    def testHt(self):
        testdemes = self.make_test_demes()
        Ht = testdemes.Ht()
        testvalues = larry.fromtuples([('Locus 1', 0.5),
                                       ('Locus 2', 0.5)])
        self.assertEqual(Ht, testvalues, 'Incorrect Ht values')
        
    def testHt_prime_est(self):
        testdemes = self.make_test_demes()
        Ht_prime_est = testdemes.Ht_prime_est()
        testvalues = larry.fromtuples([('Locus 1', 0.5), 
                                       ('Locus 2', 0.5)])
        self.assertEqual(Ht_prime_est, testvalues, "Incorrect Ht'-est values") 
        
    def testAlleleCounts(self):
        testdemes = self.make_test_demes()
        allele_counts = testdemes.allele_counts()
        testvalues = larry.fromtuples([('PopA', 'Locus 1', 10.0), 
                                       ('PopA', 'Locus 2', 10.0),
                                       ('PopB', 'Locus 1', 10.0),
                                       ('PopB', 'Locus 2', 10.0)])
        self.assertEqual(allele_counts, testvalues, 'Incorrect Allele Counts')
        
    def testLociHarmonicMeans(self):
        testdemes = self.make_test_demes()
        loci_harmonic_means = testdemes.loci_harmonic_means()
        testvalues = larry.fromtuples([('Locus 1', 10.0), 
                                       ('Locus 2', 10.0)])
        self.assertEqual(loci_harmonic_means, testvalues, 'Incorrect Harmonic Mean at Loci')
    
    def testHs_est(self):
        testdemes = self.make_test_demes()
        Hs_est = testdemes.Hs_est()
        testvalues = larry.fromtuples([('Locus 1', 0.33684210526315772),
                                       ('Locus 2', 0.52631578947368418)])
        self.assertEqual(Hs_est, testvalues, ) # 'Incorrect Hs-est values'
        
    def testHt_est(self):
        testdemes = self.make_test_demes()
        Ht_est = testdemes.Ht_est()
        testvalues = larry.fromtuples([('Locus 1', 0.508421052631579),
                                       ('Locus 2', 0.51315789473684215)])
        self.assertEqual(Ht_est, testvalues, 'Incorrect Ht-est values.')
        
    def testHst(self):
        testdemes = self.make_test_demes()
        Hst = testdemes.Hst()
        testvalues = larry.fromtuples([('Locus 1', 0.26470588235294135),
                                       ('Locus 2', 0.0)])
        self.assertEqual(Hst, testvalues, 'Incorrect Hst values.')
        
    def testDelta_s(self):
        testdemes = self.make_test_demes()
        delta_s = testdemes.delta_s()
        testvalues = larry.fromtuples([('Locus 1', 1.4705882352941173),
                                       ('Locus 2', 2.0)])
        self.assertEqual(delta_s, testvalues, 'Incorrect delata-s values')
        
    def testDelta_t(self):
        testdemes = self.make_test_demes()
        delta_t = testdemes.delta_t()
        testvalues = larry.fromtuples([('Locus 1', 2.0), 
                                       ('Locus 2', 2.0)])
        self.assertEqual(delta_t, testvalues, 'Incorrect delta-t values.')
        
    def testDst(self):
        testdemes = self.make_test_demes()
        Dst = testdemes.Dst()
        testvalues = larry.fromtuples([('Locus 1', 0.18000000000000016),
                                       ('Locus 2', 0.0)])
        self.assertEqual(Dst, testvalues, 'Incorrect Dst values.')
        
    def testGst(self):
        testdemes = self.make_test_demes()
        Gst = testdemes.Gst()
        testvalues = larry.fromtuples([('Locus 1', 0.36000000000000032),
                                       ('Locus 2', 0.0)])
        self.assertEqual(Gst, testvalues, 'Incorrect Gst values.')
        
    def testD(self):
        testdemes = self.make_test_demes()
        D = testdemes.D()
        testvalues = larry.fromtuples([('Locus 1', 0.52941176470588269),
                                       ('Locus 2', 0.0)])
        self.assertEqual(D, testvalues, 'Incorrect D values')
        
    def testGst_est(self):
        testdemes = self.make_test_demes()
        Gst_est=  testdemes.Gst_est()
        testvalues = larry.fromtuples([('Locus 1', 0.33747412008281613),
                                       ('Locus 2', -0.025641025641025501)])
        self.assertEqual(Gst_est, testvalues, 'Incorrect Gst-est values.')
        
    def testG_prime_st_est(self):
        testdemes = self.make_test_demes()
        G_prime_st_est = testdemes.G_prime_st_est()
        testvalues = larry.fromtuples([('Locus 1', 0.68030497223043851),
                                       ('Locus 2', -0.08262108262108217)])
        self.assertEqual(G_prime_st_est, testvalues, "Incorrect G'st-est values.")
        
    def testG_double_prime_st_est(self):
        testdemes = self.make_test_demes()
        G_double_prime_st_est = testdemes.G_double_prime_st_est().totuples()
        testvalues = larry.fromtuples([('Locus 1', 0.76097105508870266),
                                        ('Locus 2', -0.11111111111111048)])
        self.assertEqual(G_double_prime_st_est, testvalues, "Incorrect G''st values.")
    
    def testD_est(self):
        testdemes = self.make_test_demes()
        D_est = testdemes.D_est()
        testvalues = larry.fromtuples([('Locus 1', 0.51746031746031795),
                                       ('Locus 2', -0.055555555555555254)])
        self.assertEqual(D_est, testvalues, 'Incorrect D-est values.')
    
    def testMultilocusGst_est(self):
        testdemes = self.make_test_demes()
        multilocusGst_est = testdemes.multilocusGst_est()
        testvalue = 0.15507470376094823
        self.assertEqual(multilocusGst_est,testvalue, 'Incorrect multilocus Gst-est values.')
    
    def testMultilocusG_prime_est(self):
        testdemes = self.make_test_demes()
        multilocusG_prime_est = testdemes.multilocusG_prime_st_est()
        testvalue = 0.39055851317572143
        self.assertEqual(multilocusG_prime_est, testvalue, "Incorrect multilocus G'st value.")
    
    def testMultilocusG_double_prime_st_est(self):
        testdemes = self.make_test_demes()
        multilocusG_double_prime_st_est = testdemes.multilocusG_double_prime_st_est()
        testvalue =  0.47237915881983739
        self.assertEqual(multilocusG_double_prime_st_est, testvalue, "Incorrect multilocus G''st-est value.")
        
    def testMultilocusD_est(self):
        testdemes = self.make_test_demes()
        multilocusD_est = testdemes.multilocusD_est()
        testvalue = 0.090963240209061477
        self.assertEqual(multilocusD_est, testvalue, 'Incorrect multilocus D-est value.')

class MonteCarloTests(unittest.TestCase):

        #        Locus 1    Locus 2
    popA = [['001001', '002002'],
            ['001001', '002002'],
            ['002002', '002002'],
            ['002002', '002002'],
            ['002002', '002002'],
            ['002002', '001001'],
            ['002002', '001001'],
            ['002002', '001001'],
            ['002002', '001001'],
            ['002002', '001001']]

    popB = [['002002', '003003'],
            ['002002', '003003'],
            ['001001', '003003'],
            ['001001', '003003'],
            ['001001', '003003'],
            ['001001', '001001'],
            ['001001', '001001'],
            ['001001', '001001'],
            ['001001', '001001'],
            ['001001', '001001']]
            
    def make_test_demes(self):
        """setups up basic test deme"""
        pop1 = population(self.popA, loci=['Locus 1', 'Locus 2'], name = 'PopA')
        pop2 = population(self.popB, loci=['Locus 1', 'Locus 2'], name = 'PopB')
        testdemes = populations()
        testdemes.append(pop1)
        testdemes.append(pop2)
        return testdemes


    def testSimulate(self):
        demes = self.make_test_demes()
        MtC = MonteCarlo()
        sim_data = MtC.simulate(demes, estimator="D_est", replicates=50)
        MtC.p_values(sim_data, demes, estimator="D_est")
        pass


if __name__ == '__main__':
    
    args = get_args()
    fin = open(args.input_file, 'rU')
    lines = fin.readlines()
    demes = parse_genepop(lines)
    fin.close()
    
    pops = parse_genepop(lines)
    MtC = MonteCarlo()
    estimator = "G_prime_st_est"
    sim_data = MtC.simulate(demes, estimator=estimator, replicates=50)
    MtC.p_values(sim_data, demes, estimator=estimator)
    
    
    


    # new_pair = populations()
    # new_pair.append(pops.pops[4])
    # new_pair.append(pops.pops[-1])
    # print new_pair.pop_names()
    # print new_pair.multilocusGst_est()
    # pops.pairwise_estimators('Gst_est')
    # pops.pairwise_estimators('G_prime_st_est')
    # pops.pairwise_estimators('G_double_prime_st_est')
    # pops.pairwise_estimators('D_est')
    
    
    #data =  pops.larry2csv()
    # #print 'Ht', pops.Ht()
    # # print 'loci_harmonic_means', pops.loci_harmonic_means()
    #test = pops.__allele_freqs_3D_larry__()
    # alleles = pops.pops[0].allele_freqs()
    # print 'Hs', pops.Hs_est()
    # print 'Hs_prime_est', pops.Hs_prime_est()
    # print 'Ht_prime_est', pops.Ht_prime_est()
    # print pops.all_estimators()
    # pops.all_multilocus_estimators()
    # # print 'Gst_est', pops.Gst_est().totuples()
    #print 'Gst_est', pops.multilocusGst_est()
    # # print "G'st_est", pops.G_prime_st_est().totuples()
    # # print "G'st_est", pops.multilocusG_prime_st_est()
    # print "G''st_est", pops.G_double_prime_st_est().totuples()
    # print "G''st_est", pops.multilocusG_double_prime_st_est()
    # print "D_est", pops.D_est().totuples()
    # print "D_est", pops.multilocusD_est()
    #result = pops.generate_null_distribution(estimator='Gst_est', replicates=5)
    
    #obs = pops.calculate_pairwise_estimators(estimator='Gst_est')

    # suite = unittest.TestLoader().loadTestsFromTestCase(InputFileTest)
    # unittest.TextTestRunner(verbosity=5).run(suite)
    
    # suite = unittest.TestLoader().loadTestsFromTestCase(PopulationTests)
    # unittest.TextTestRunner(verbosity=1).run(suite)
    # suite = unittest.TestLoader().loadTestsFromTestCase(PopulationsTests)
    # unittest.TextTestRunner(verbosity=1).run(suite)

    # suite = unittest.TestLoader().loadTestsFromTestCase(MonteCarloTests)
    # unittest.TextTestRunner(verbosity=1).run(suite)
