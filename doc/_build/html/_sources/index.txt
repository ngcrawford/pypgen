.. pypgen documentation master file, created by
   sphinx-quickstart on Sun Feb  3 11:35:09 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to pypgen's documentation!
==================================

Pypgen provides various utilities for estimating standard genetic diversity measures including Gst,
G'st, G''st, and Jost's D from large genomic datasets (Hedrick, 2005; Jost, 2008; Masatoshi Nei,
1973; Nei & Chesser, 1983). Pypgen operates both on individual SNPs as well as on user defined
regions (e.g., five kilobase windows tiled across each chromosome). For the windowed analyses
pypgen estimates the multi-locus versions of each estimator.


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


Prerequisites:
++++++++++++++

Pypgen is written in Python 2.7. It may run under Python 2.6, but I haven't tested it. It doesn't run under Python 3. In order to interact with bgziped files it requires `samtools <http://samtools.sourceforge.net/>`_ and `pysam <http://www.cgat.org/~andreas/documentation/pysam/contents.html>`_ to be installed.


Quick Installation:
+++++++++++++++++++

If you already have a working install of pysam, pypgen can be installed from `PyPi <http://pypi.python.org/pypi/pypgen>`_  using `pip <http://pypi.python.org/pypi/pip>`_ or `setuptools <http://pypi.python.org/pypi/setuptools>`_:

::

        pip install pypgen

or, 

::

        easy_install -U pypgen

However, it's recommended, at least in these early days of pypgen, to install it directly from the github repository:


::

       pip install -e git+https://github.com/ngcrawford/pypgen.git#egg=Package
	   
 
Contents:
+++++++++

	.. toctree::
   	   :maxdepth: 2

   	   install


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

