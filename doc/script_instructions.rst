Scripts:
========

vcf_snpwise_fstats
++++++++++++++++++

This script calculates *F*-statistics for each pair of populations at each SNV in the supplied region.  

**Working Example:**

    Note that ``path/to/pypgen/data/example.vcf.gz`` needs to be updated to the directory in which the source code for ``pypgen`` is found.


    ::
    
        vcf_snpwise_fstats.py \
          -i pypgen/data/example.vcf.gz \
          -p outgroups:h665,i02-210 \
             pop1:c511,c512,c513,c514,c515,c563,c614,c630,c639,c640 \
             pop2:m523,m524,m525,m589,m675,m676,m682,m683,m687,m689 \
          -c 2 \
          -r Chr01:1-10001 | head


Command Line Flags  
------------------


**Input:** [ ``-i``, ``--input`` ]

    Defines the path to the input VCF file.

**Output:** [ ``-o``, ``--output`` ]

    Defines the path to the output csv/txt file. If it's not set it defaults to standard out (stout).

**Cores:** [ ``-c``, ``--cores`` ]

    The number of cores to use.

**Regions:** [ ``-r``, ``-R``, ``--regions`` ]

    This allows for selecting a subset of the VCF file for analysis. The command format should familiar to if you use GATK or samtools. A region can be presented, for example, in the following format: ‘chr2’ (the whole chr2), ‘chr2:1000000’ (region starting from 1,000,000bp) or ‘chr2:1,000,000-2,000,000’ (region between 1,000,000 and 2,000,000bp including the end points). The coordinate system is 1-based. Multiple regions can be submitted separated by spaces. [Note: this is the same format as samtools/GATK and this example text is largely borrowed from samtools]


**Populations:** [ ``-p``, ``--populations`` ]

    Names of populations and samples. The format is: "PopName:sample1,sample2,.. PopName2:sample3,sample4,..." with colons after each population name and samples delimited by commas. Whitespace is used to delimit populations.
    

**Minimum Number of Samples:** [ ``-m``, ``--min-samples`` ]

    This allows one to set the minimum number of samples per population that a SNV needs to have in order to be included in the analysis.
    

**Column Separator:** [ ``-s``, ``--column-separator`` ]

    This allows one to set the separator to be uses in the output. The default value is ``,`` which makes the output comma separated (csv). If you're planning on using tabix to index the output you'll need to set the sep to ``\t``.
              
**Zero Based:** [ ``--zero-based`` ]

    Setting this flag makes the output positions zero based (e.g., BED like).   


Output 
------


- The chrom and pos columns are fixed in positions 1 and 2, but the rest of the columns vary depending on the number of populations being compared and their names.


    +---------------------------------------+-------------------------------------------------+
    | Label:                                | Definition:                                     |
    +=======================================+=================================================+
    | *chrom*                               | ID of chromosome/scaffold/contig/etc.           |
    +---------------------------------------+-------------------------------------------------+
    | *pos*                                 | Position of SNP                                 |
    +---------------------------------------+-------------------------------------------------+
    | *pop1.sample_count*                   | Number of samples represented                   |
    +---------------------------------------+-------------------------------------------------+
    | cont.                                 | Additional population sample counts             |
    +---------------------------------------+-------------------------------------------------+
    | *Pop1.Pop2.D\_est*\                   | D corrected for sample size (Jost 2008)         |
    +---------------------------------------+-------------------------------------------------+
    | *Pop1.Pop2.D\_est.stdev*\             | Standard deviation D corrected for sample size  |
    +---------------------------------------+-------------------------------------------------+
    | *Pop1.Pop2.G\_double\_prime\_st\_est* | Corrected Hedrick’s G'st                        |
    |                                       | (Meirmans & Hedrick 2011)                       |
    +---------------------------------------+-------------------------------------------------+
    | *Pop1.Pop2.G\_prime\_st\_est*         | Standardized Gst (Hedrick 2005)                 |
    +---------------------------------------+-------------------------------------------------+
    | *Pop1.Pop2.Gst\_est*                  | Fst corrected for sample size and allowing for  |
    |                                       | multiallelic loci (Nei & Chesser 1983)          |
    +---------------------------------------+-------------------------------------------------+
    | *Pop1.Pop2.Hs\_est*                   | Within-population gene/locus diversity          |
    |                                       | (e.g., expected heterozygosity)                 |
    +---------------------------------------+-------------------------------------------------+
    | *Pop1.Pop2.Ht\_est*                   | Total gene/locus diversity                      |
    +---------------------------------------+-------------------------------------------------+
    | cont...                               | Pairwise comparisons of F-statistics cont...    |
    +---------------------------------------+-------------------------------------------------+
    |*Pop1\_fixed*                          | If a sample is fixed at a particular allele     |
    |                                       | this flag is set to 1 (= "True" in binary)      |    
    +---------------------------------------+-------------------------------------------------+
    | cont...                               | Additional fixed SNPs cont...                   |
    +---------------------------------------+-------------------------------------------------+



vcf_sliding_window
++++++++++++++++++

This script calculates *F*-statistics for each pair of populations at each window in the supplied region. This script requires that the input VCF file be bgzipped because it uses ``tabix`` to extract the windows. 

**Working Example:**

    Note that ``path/to/pypgen/data/example.vcf.gz`` needs to be updated to the directory in which the source code for ``pypgen`` is found.

    ::

        vcf_sliding_window.py \
          -i path/to/pypgen/data/example.vcf.gz \
          -p outgroups:h665,i02-210 \
             pop1:c511,c512,c513,c514,c515,c563,c614,c630,c639,c640 \
             pop2:m523,m524,m525,m589,m675,m676,m682,m683,m687,m689 \
          -c 2 \
          -w 5000 \
          -r Chr01:1-10001 | head


Command Line Flags  
------------------

``vcf_sliding_window.py`` shares the same commands as ``vcf_snpwise_fstats.py`` with the single addition of a window size flag.


**Window Size:** [ ``-w``, ``--window-size`` ]

    Windows are non overlapping and start at the first bp in the particular chromosome. 



Output 
------

**vcf\_sliding\_window.py:** 

- The format is loosely based on the `BED specification <http://genome.ucsc.edu/FAQ/FAQformat.html#format1>`_. Although the first three column IDs will remain static for the foreseeable future, I expect to add more fields as I add additional functionality to pypgen. Also, the default output is one based, but it is possible to make the positions zero based by including the ``--zero-based`` flag when you run the script.

- The population IDs and the total number of populations come from those defined by the user. This means the number of pairwise population comparisons and hence the total number of columns is conditional on the number of defined populations. 

    +---------------------------------------------+-----------------------------------------------+
    | Label:                                      | Definition:                                   |
    +=============================================+===============================================+
    | *chrom*                                     | ID of chromosome/scaffold/contig/etc.         |
    +---------------------------------------------+-----------------------------------------------+
    | *chromStart*                                | Starting position of window                   |
    +---------------------------------------------+-----------------------------------------------+
    | *chromEnd*                                  | Ending position of window                     |
    +---------------------------------------------+-----------------------------------------------+
    | *snp\_count*                                | Total Number of SNPs in window                |
    +---------------------------------------------+-----------------------------------------------+
    | *total\_depth\_mean*                        | Mean read depth across window                 |
    +---------------------------------------------+-----------------------------------------------+
    | *total\_depth\_stdev*                       | Standard deviation of read depth across window|
    +---------------------------------------------+-----------------------------------------------+
    | *Pop1.sample\_count.mean*                   | Mean number of samples per snp for 'Pop1'     |
    +---------------------------------------------+-----------------------------------------------+
    | *Pop1.sample\_count.stdev*                  | Standard deviation of samples per snp for     |
    |                                             | 'Pop1'                                        |
    +---------------------------------------------+-----------------------------------------------+
    | *Pop2.sample\_count.mean*                   | Mean number of samples per snp for 'Pop2'     |
    +---------------------------------------------+-----------------------------------------------+
    | *Pop2.sample\_count.stdev*                  | Standard deviation of samples per snp for     |
    |                                             | 'Pop2'                                        |
    +---------------------------------------------+-----------------------------------------------+
    | *Pop2.Pop1.D\_est*                          | Multilocus D_est (Jost 2008)                  |
    +---------------------------------------------+-----------------------------------------------+
    | *Pop2.Pop1.D\_est.stdev*                    | Standard Deviation of SNVwise D_est across    |
    |                                             | the window                                    |
    +---------------------------------------------+-----------------------------------------------+
    | *Pop2.Pop1.G\_double\_prime\_st\_est*       | Corrected Hedrick’s G'st                      |
    |                                             | (Meirmans & Hedrick 2011)                     |
    +---------------------------------------------+-----------------------------------------------+
    | *Pop2.Pop1.G\_double\_prime\_st\_est.stdev* | Standard Deviation of Corrected Hedrick’s     | 
    |                                             | G'st across the window                        |
    +---------------------------------------------+-----------------------------------------------+
    | *Pop2.Pop1.G\_prime\_st\_est*               | Standardized Gst (Hedrick 2005)               |
    +---------------------------------------------+-----------------------------------------------+
    | *Pop2.Pop1.G\_prime\_st\_est.stdev*         | Standard Deviation of Standardized Gst across |
    |                                             | the window                                    |
    +---------------------------------------------+-----------------------------------------------+
    | *Pop2.Pop1.Gst\_est*                        | Fst corrected for sample size and             |
    |                                             | allowing for multiallelic loci                |
    |                                             | (Nei & Chesser 1983)                          |
    +---------------------------------------------+-----------------------------------------------+
    | *Pop2.Pop1.Gst\_est.stdev*                  | Standard Deviation of Fst corrected for       |
    |                                             | sample size and  allowing for multiallelic    |
    |                                             | loci (Nei & Chesser 1983)                     |
    +---------------------------------------------+-----------------------------------------------+
    | cont...                                     | The rest of the pairwise comparisons follow...|
    +---------------------------------------------+-----------------------------------------------+


