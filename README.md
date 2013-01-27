Welcome to Pypgen (v0.2.0) *BETA*
--------------------------

Pypgen provides various utilities for estimating standard genetic diversity measures including Gst, G'st, G''st, and Jost's D from large genomic datasets (Hedrick, 2005; Jost, 2008; Masatoshi Nei, 1973; Nei & Chesser, 1983). Pypgen operates both on the level of individual SNPs as well as on user defined regions (e.g., five kilobase windows tiled across each chromosome). For the windowed analyses, pypgen estimates the multi-locus versions of each estimators.

### Features:
- Operates on standard [VCF (Variant Call Format)][1] formatted SNP calls
- Uses [bgziped][5] input for fast random access
- Takes advantage of multiple processor cores
- Handles multiallelic SNP calls
- Allows a single VCF file to contain multiple populations
- Calculates additional metrics:
	- snp count
	- populations with fixed alleles (per SNP)
    - mean read depth (+/- STDEV) 
    - more as I think of them

### Important Note:
PYPGEN IS STILL IN ACTIVE DEVELOPMENT AND ALMOST CERTAINLY CONTAINS BUGS. 
If you find a bug please file a report in the 'issues' section of this repository and I'll address it as soon as I can. 

### Enclosed Scripts:
- Sliding window analysis (vcf_sliding_window.py) 
- Per SNP analysis (vcf_snpwise_fstats.py)
- UnitTests (vcf_tests.py)

### Dependancies:
- OSX or Linux
- [Python 2.7][4]
- [Numpy][2]
- [pysam][3] and [samtools][8]

Note: a *setup.py* script that will automate the installation of the depencancies will be forthcoming. 

### Installation:
First install the dependancies. I like to use [pip][6] for this purpose. If running on OS X I recommend using [homebrew][7] to install samtools, requirement of pysam, prior to installing pysam. Once the dependancies are installed clone this repository:
        
        git clone git@github.com:ngcrawford/pypgen.git

You should be able to run the UnitTests without any problems:

        python vcf_tests.py

Infomation about each script can be obtained by running:

        python [script name].py -h

### Citations:
- Crawford, N. G. (2009). SMOGD: software for the measurement of genetic diversity. Molecular Ecology Resources, 10, 556-557. doi:10.1111/j.1755-0998.2009.02801.x
- Hedrick, P. W. (2005). A standardized genetic differentiation measure. Evolution; international journal of organic evolution, 59(8), 1633-1638.
- Jost, L. (2008). Gst and its relatives do not measure differentiation. Molecular Ecology, 17(18), 4015-4026. doi:10.1111/mec.2008.17.issue-18
- Nei, M. (1973). Analysis of Gene Diversity in Subdivided Populations. Proceedings of the National Academy of Sciences of the United States of America, 70(12 Pt 1-2), 3321.
- Nei, M., & Chesser, R. K. (1983). Estimation of fixation indices and gene diversities. Annals of human genetics, 47(Pt 3), 253-259.


[1]: http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41
[2]: http://www.numpy.org
[3]: http://wwwfgu.anat.ox.ac.uk/~andreas/documentation/samtools/contents.html
[4]: http://www.python.org/download/releases/2.7/
[5]: http://samtools.sourceforge.net/tabix.shtml
[6]: http://pypi.python.org/pypi/pip
[7]: http://mxcl.github.com/homebrew/
[8]: http://samtools.sourceforge.net/
