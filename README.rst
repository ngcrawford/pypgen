Welcome to Pypgen (v0.2.0) *BETA*
---------------------------------

Pypgen provides various utilities for estimating standard genetic
diversity measures including Gst, G'st, G''st, and Jost's D from large
genomic datasets (Hedrick, 2005; Jost, 2008; Masatoshi Nei, 1973; Nei &
Chesser, 1983). Pypgen operates both on the level of individual SNPs as
well as on user defined regions (e.g., five kilobase windows tiled
across each chromosome). For the windowed analyses, pypgen estimates the
multi-locus versions of each estimators.

Features:
~~~~~~~~~

-  Operates on standard `VCF (Variant Call
   Format) <http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41>`_
   formatted SNP calls
-  Uses `bgziped <http://samtools.sourceforge.net/tabix.shtml>`_ input
   for fast random access
-  Takes advantage of multiple processor cores
-  Handles multiallelic SNP calls
-  Allows a single VCF file to contain multiple populations
-  Calculates additional metrics:

   -  snp count per window
   -  mean read depth (+/- STDEV) per window

   -  populations with fixed alleles per SNP
   -  more as I think of them

Important Note:
~~~~~~~~~~~~~~~

PYPGEN IS STILL IN ACTIVE DEVELOPMENT AND ALMOST CERTAINLY CONTAINS
BUGS. If you find a bug please file a report in the 'issues' section of
this repository and I'll address it as soon as I can.

Enclosed Scripts:
~~~~~~~~~~~~~~~~~

-  Sliding window analysis (vcf\_sliding\_window.py)
-  Per SNP analysis (vcf\_snpwise\_fstats.py)

Dependancies:
~~~~~~~~~~~~~

-  OSX or Linux
-  `Python 2.7 <http://www.python.org/download/releases/2.7/>`_
-  `Numpy <http://www.numpy.org>`_
-  `pysam <http://wwwfgu.anat.ox.ac.uk/~andreas/documentation/samtools/contents.html>`_
   and `samtools <http://samtools.sourceforge.net/>`_

Installation:
~~~~~~~~~~~~~

First install the dependancies. I like to use
`pip <http://pypi.python.org/pypi/pip>`_ for this purpose. If running on
OS X I recommend using `homebrew <http://mxcl.github.com/homebrew/>`_ to
install samtools, requirement of pysam, prior to installing pysam. Once
the dependancies are installed clone this repository:

::

        git clone git@github.com:ngcrawford/pypgen.git

You should be able to run the UnitTests without any problems:

::

        python scr/tests/tests.py

or, install and run nose

::

        pip install nose
        nosetests

Infomation about each script can be obtained by running:

::

        python [script name].py -h

Output: (note this will probably change)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**vcf\_sliding\_window.py:**

-  *chrm* = Name of chromosome
-  *start* = Starting position of window
-  *stop* = Ending position of window
-  *snp\_count* = Total Number of SNPs in window
-  *total\_depth\_mean* = Mean read depth across window
-  *total\_depth\_stdev* = Standard deviation of read depth across
   window
-  *Pop1.sample\_count.mean* = Mean number of samples per snp for 'Pop1'
-  *Pop1.sample\_count.stdev* = Standard deviation of samples per snp
   for - 'Pop1'
-  *Pop2.sample\_count.mean* = Mean number of samples per snp for 'Pop2'
-  *Pop2.sample\_count.stdev* = Standard deviation of samples per snp
   for 'Pop2'
-  *Pop2.Pop1.D\_est* = Multilocus Dest (Jost 2008)
-  *Pop2.Pop1.G\_double\_prime\_st\_est* = (Meirmans & Hedrick 2011)
-  *Pop2.Pop1.G\_prime\_st\_est* = Standardized Gst (Hedrick 2005)
-  *Pop2.Pop1.Gst\_est* = Fst corrected for sample size and allowing for
   multiallelic loci (Nei & Chesser 1983)
-  cont...

**vcf\_snpwise\_fstats.py:**

-  *chrm* = Name of chromosome
-  *pos* = Position of SNP
-  *outgroups* = Number of samples
-  *pop1* = Population ID
-  *pop1.outgroups.D\_est*\ = Multilocus Dest (Jost 2008)
-  *pop1.outgroups.G\_double\_prime\_st\_est* = (Meirmans & Hedrick
   2011)
-  *pop1.outgroups.G\_prime\_st\_est* = Standardized Gst (Hedrick 2005)
-  *pop1.outgroups.Gst\_est* = Fst corrected for sample size and
   allowing for multiallelic loci (Nei & Chesser 1983)
-  *pop1.outgroups.Hs\_est*
-  *pop1.outgroups.Ht\_est*
-  cont...,
-  *outgroups\_fixed* = If a sample is fixed at a particular allele this
   flag is set to 1 (= "True" in binary).
-  cont...

Citations:
~~~~~~~~~~

-  Crawford, N. G. (2009). SMOGD: software for the measurement of
   genetic diversity. Molecular Ecology Resources, 10, 556-557.
   doi:10.1111/j.1755-0998.2009.02801.x
-  Hedrick, P. W. (2005). A standardized genetic differentiation
   measure. Evolution; international journal of organic evolution,
   59(8), 1633-1638.
-  Jost, L. (2008). Gst and its relatives do not measure
   differentiation. Molecular Ecology, 17(18), 4015-4026.
   doi:10.1111/mec.2008.17.issue-18
-  Meirmans, P., & Hedrick, P. (2011). Assessing population structure:
   FST and related measures. Molecular Ecology, 11, 5–8.
-  Nei, M. (1973). Analysis of Gene Diversity in Subdivided Populations.
   Proceedings of the National Academy of Sciences of the United States
   of America, 70(12 Pt 1-2), 3321.
-  Nei, M., & Chesser, R. K. (1983). Estimation of fixation indices and
   gene diversities. Annals of human genetics, 47(Pt 3), 253-259.

