Welcome to Pypgen (v0.2.0)
--------------------------

Pypgen provides various utilities for estimating standard genetic diversity measures including Gst, G'st, G''st, Jost's D from large genomic datasets (Hedrick, 2005; Jost, 2008; Masatoshi Nei, 1973; Nei & Chesser, 1983). Pypgen operates both on the level of individual SNPs as well as on user defined windows (e.g., 5 kilobases tiled across each chromosome). For the windowed analyses, pypgen estimates the multilocus versions of each of the estimators. Pypgen operates on standard [VCF (Variant Call Format)][1] formated SNP calls

Scripts:

- Sliding window analysis (vcf_sliding_window.py) 
- Per SNP analysis (vcf_snpwise_fstats.py)


Citations:

- Crawford, N. G. (2009). SMOGD: software for the measurement of genetic diversity. Molecular Ecology Resources, 10, 556-557. doi:10.1111/j.1755-0998.2009.02801.x
- Hedrick, P. W. (2005). A standardized genetic differentiation measure. Evolution; international journal of organic evolution, 59(8), 1633-1638.
- Jost, L. (2008). Gst and its relatives do not measure differentiation. Molecular Ecology, 17(18), 4015-4026. doi:10.1111/mec.2008.17.issue-18
- Nei, M. (1973). Analysis of Gene Diversity in Subdivided Populations. Proceedings of the National Academy of Sciences of the United States of America, 70(12 Pt 1-2), 3321.
- Nei, M., & Chesser, R. K. (1983). Estimation of fixation indices and gene diversities. Annals of human genetics, 47(Pt 3), 253-259.

[1]: http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41