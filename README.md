Welcome to Pypgen (v0.2.0)
--------------------------

Pypgen provides various utilities for estimating standard genetic diversity measures such as Gst, G'st, G''st, and Jost's D from large genomic datasets (Hedrick, 2005; Jost, 2008; Masatoshi Nei, 1973; Nei & Chesser, 1983). Pypgen operates both on the level of individual SNPs as well as on user defined regions (e.g., five kilobase windows tiled across each chromosome). For the windowed analyses, pypgen estimates the multi-locus versions of each of the estimators. Pypgen operates on standard [VCF (Variant Call Format)][1] formatted SNP calls.

### Features:
- Takes advantage of multiple processor cores
- Handles multiallelic SNP calls
- Allows a single VCF file to contain multiple populations
- Calculates additional summary statistics:
	- snp count 
	- mean read depth and standard deviation 
- UnitTests

### Important Notes:
PYPGEN IS STILL IN ACTIVE DEVELOPMENT AND ALMOST CERTAINLY CONTAINS BUGS. IF YOU FIND A BUG PLEASE FILE A REPORT IN THE 'ISSUES' SECTION AND I'LL ADDRESS IT AS SOON AS I CAN. 

### Enclosed Scripts:

- Sliding window analysis (vcf_sliding_window.py) 
- Per SNP analysis (vcf_snpwise_fstats.py)
- UnitTests (vcf_tests.py)

### Dependancies:

- Python 2.7
- Numpy 

### Citations:
<small>
- Crawford, N. G. (2009). SMOGD: software for the measurement of genetic diversity. Molecular Ecology Resources, 10, 556-557. doi:10.1111/j.1755-0998.2009.02801.x
- Hedrick, P. W. (2005). A standardized genetic differentiation measure. Evolution; international journal of organic evolution, 59(8), 1633-1638.
- Jost, L. (2008). Gst and its relatives do not measure differentiation. Molecular Ecology, 17(18), 4015-4026. doi:10.1111/mec.2008.17.issue-18
- Nei, M. (1973). Analysis of Gene Diversity in Subdivided Populations. Proceedings of the National Academy of Sciences of the United States of America, 70(12 Pt 1-2), 3321.
- Nei, M., & Chesser, R. K. (1983). Estimation of fixation indices and gene diversities. Annals of human genetics, 47(Pt 3), 253-259.
</small>

[1]: http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41