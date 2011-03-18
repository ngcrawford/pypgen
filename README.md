Welcome to Pypgen (v0.1)
-----------------------

This is repository contains a python classes for parsing genepop files calculating *Gst*, *G'st*, *G''st*, *D*, and the 'estimators' of each of these diversity measures (Hedrick, 2005; Jost, 2008; Masatoshi Nei, 1973; Nei & Chesser, 1983). This code is a rewrite of my [SMOGD](https://github.com/ngcrawford/SMOGD) code.  Although, at this time, there are only minimal changes to the functionality, the code is considerably more pythonic with unittests and classes.  This should make it substantially easier to maintain. I'm also hoping that population geneticists will fork this code and add more functionality. 

### Up Next:###

1.) **Faster calculations:** Currently each diversity measure is calculated from the raw genotypes each time its calculator function is called. This is slow if you wish to calculate more than one diversity measure at a time. Allowing function calls to use pre-calculated data should speed thing up considerably. For example, calculating *Gst* requires calculating *Hs* and *Ht*.  The way the code is written now calculating a subsequent calculation *D* involves re-calculating *Hs* and *Ht*. In the future I'll add code to store the relevant intermediate steps.

2.) **Missing Functionality:** The pairwise calculators, bootstrapping analysis, and Arlequin format parser didn't make this version.

3.) **New Functionality:** One of the main goals of this rewrite is to provide command-line tools.  I'll be producing a command-line version of shortly. I'm also planning on including the *D-est-chao* estimator (eq 13 from Jost, 2008).


Hedrick, P. W. (2005). A standardized genetic differentiation measure. Evolution; international journal of organic evolution, 59(8), 1633-1638.
Jost, L. (2008). Gst and its relatives do not measure differentiation. Molecular Ecology, 17(18), 4015-4026. doi:10.1111/mec.2008.17.issue-18
Masatoshi Nei. (1973). Analysis of Gene Diversity in Subdivided Populations. Proceedings of the National Academy of Sciences of the United States of America, 70(12 Pt 1-2), 3321.
Nei, M., & Chesser, R. K. (1983). Estimation of fixation indices and gene diversities. Annals of human genetics, 47(Pt 3), 253-259.
