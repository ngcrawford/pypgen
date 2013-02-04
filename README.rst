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
-  `pysam <http://wwwfgu.anat.ox.ac.uk/+andreas/documentation/samtools/contents.html>`_
   and `samtools <http://samtools.sourceforge.net/>`_


Documentation:
++++++++++++++

Detailed documentation is available on `ReadTheDocs <https://pypgen.readthedocs.org/en/latest/index.html>`_. It includes a tutorial and installation instructions.

