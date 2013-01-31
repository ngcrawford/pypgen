Welcome to Pypgen (v0.2.1) *BETA*
---------------------------------

Pypgen provides various utilities for estimating standard genetic
diversity measures including Gst, G'st, G''st, and Jost's D from large
genomic datasets (Hedrick, 2005; Jost, 2008; Masatoshi Nei, 1973; Nei &
Chesser, 1983). Pypgen operates both on individual SNPs as
well as on user defined regions (e.g., five kilobase windows tiled
across each chromosome). For the windowed analyses pypgen estimates the
multi-locus versions of each estimator.

Features:
+++++++++

-  Handles multiallelic SNP calls
-  Allows a single VCF file to contain multiple populations
-  Operates on standard `VCF (Variant Call
   Format) <http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41>`_
   formatted SNP calls
-  Uses `bgziped <http://samtools.sourceforge.net/tabix.shtml>`_ input
   for fast random access
-  Takes advantage of multiple processor cores
-  Calculates additional metrics:

   -  snp count per window
   -  mean read depth (+/- STDEV) per window
   -  populations with fixed alleles per SNP
   -  more as I think of them

Important Note:
+++++++++++++++

PYPGEN IS STILL IN ACTIVE DEVELOPMENT AND ALMOST CERTAINLY CONTAINS
BUGS. If you find a bug please file a report in the `issues section <https://github.com/ngcrawford/pypgen/issues>`_ of
the github repository and I'll address it as soon as I can.

Enclosed Scripts:
+++++++++++++++++

-  Sliding window analysis (vcf\_sliding\_window.py)
-  Per SNP analysis (vcf\_snpwise\_fstats.py)

Dependancies:
+++++++++++++

-  OSX or Linux
-  `Python 2.7 <http://www.python.org/download/releases/2.7/>`_
-  `Numpy <http://www.numpy.org>`_
-  `pysam <http://wwwfgu.anat.ox.ac.uk/+andreas/documentation/samtools/contents.html>`_
   and `samtools <http://samtools.sourceforge.net/>`_

Installation:
+++++++++++++

First install `samtools <http://samtools.sourceforge.net/>`_. On OS X I recommend using `homebrew <http://mxcl.github.com/homebrew/>`_ to do this. Once you have samtools installed and available in terminal you can use either pip or setuptools to install the current release of pypgen:

::

        pip install pypgen

or, 

::

        easy_install pypgen


of course it is still possible to install directly from the source via:

::

		tar -xzvf ~/your/download/location/pypgen-*.tar.gz
		python setup.py install

Alternately, if you like to live on the edge, you can clone and install the current development version from github.

::

       pip install -e git+https://github.com/ngcrawford/pypgen.git#egg=Package

Documentation:
++++++++++++++

More detailed documentation will be forthcoming, but in the meantime information about each script can be obtained by running:

::

        [script name].py --help 

Output: 
+++++++
**vcf\_sliding\_window.py:** 

- The format is loosely based on the `BED specification <http://genome.ucsc.edu/FAQ/FAQformat.html#format1>`_. Although the first three column IDs will remain static for the foreseeable future, I expect to add more fields as I add additional functionality to pypgen. 

- Currently positional values are one based, but this may change as this is not 'BED-like' (not zero based). 

- The population IDs and the total number of populations come from those defined by the user. This means the number of pairwise population comparisons and hence the total number of columns is conditional on the number of defined populations. 

+---------------------------------------+-------------------------------------------------+
| Label:                                | Definition:                                     |
+=======================================+=================================================+
| *chrom*                               | ID of chromosome/scaffold/contig/etc.           |
+---------------------------------------+-------------------------------------------------+
| *chromStart*                          | Starting position of window                     |
+---------------------------------------+-------------------------------------------------+
| *chromEnd*                            | Ending position of window                       |
+---------------------------------------+-------------------------------------------------+
| *snp\_count*                          | Total Number of SNPs in window                  |
+---------------------------------------+-------------------------------------------------+
| *total\_depth\_mean*                  | Mean read depth across window                   |
+---------------------------------------+-------------------------------------------------+
| *total\_depth\_stdev*                 | Standard deviation of read depth across window  |
+---------------------------------------+-------------------------------------------------+
| *Pop1.sample\_count.mean*             | Mean number of samples per snp for 'Pop1'       |
+---------------------------------------+-------------------------------------------------+
| *Pop1.sample\_count.stdev*            | Standard deviation of samples per snp for 'Pop1'|
+---------------------------------------+-------------------------------------------------+
| *Pop2.sample\_count.mean*             | Mean number of samples per snp for 'Pop2'       |
+---------------------------------------+-------------------------------------------------+
| *Pop2.sample\_count.stdev*            | Standard deviation of samples per snp for 'Pop2'|
+---------------------------------------+-------------------------------------------------+
| *Pop2.Pop1.D\_est*                    | Multilocus D (Jost 2008)                        |
+---------------------------------------+-------------------------------------------------+
| *Pop2.Pop1.G\_double\_prime\_st\_est* | Corrected Hedrick’s G'st                        |
|                                       | (Meirmans & Hedrick 2011)                       |
+---------------------------------------+-------------------------------------------------+
| *Pop2.Pop1.G\_prime\_st\_est*         | Standardized Gst (Hedrick 2005)                 |
+---------------------------------------+-------------------------------------------------+
| *Pop2.Pop1.Gst\_est*                  | Fst corrected for sample size and               |
|                                       | allowing for multiallelic loci                  |
|                                       | (Nei & Chesser 1983)                            |
+---------------------------------------+-------------------------------------------------+
| cont...                               | The rest of the pairwise comparisons follow...  |
+---------------------------------------+-------------------------------------------------+

**vcf\_snpwise\_fstats.py:**

- The chrom and pos column are fixed in positions 1 and 2, but the rest of the columns are 


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