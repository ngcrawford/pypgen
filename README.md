Welcome to Pypgen (v0.1)
-----------------------

This is repository contains a python classes for parsing [genepop files](http://genepop.curtin.edu.au/help_input.html), calculating *Gst*, *G'st*, *G''st*, *D*, and calculating the 'estimators' and multilocus parameters of each of these diversity measures (Hedrick, 2005; Jost, 2008; Masatoshi Nei, 1973; Nei & Chesser, 1983). This code is a rewrite of my [SMOGD](https://github.com/ngcrawford/SMOGD) code (Crawford, 2009) although, at this time, there are only minimal changes to the functionality.  Rather pypgen is considerably more pythonic and object oriented with unittests and classes. This should make it substantially easier to maintain. I'm also hoping that other python coding population geneticists will fork this code and add more functionality. 

### Installation: ###

Pypop requires that you install the python mathematics and plotting software [MatPlotLib](http://matplotlib.sourceforge.net/index.html).  If you're installing on OSX or linux I recommend using [easy install](http://pypi.python.org/pypi/setuptools) to do the installation as MatPlotLib can be difficult to install from source.  With that in mind the basic steps would be as follows:

1.) Download and install easy install (aka setuptools). Instructions [here](http://pypi.python.org/pypi/setuptools).

2.) Open terminal and type the following command: `easy_install matplotlib`. This should install all the dependancies for pypgen.

3.) Check that everything worked by running `python` at the command line and once in the python shell typing `import pylab`.  If you don't see any errors you're good to go. 

3.) If you're just going to run the provided command-line scripts you can jump right in, but if you're going to 'play with the code' I recommend installing ipython with the command `easy_install ipython`. You can then run a nice interactive version of python and, if you type `ipython -pylab wx` at the prompt, it'll run preloaded with all the pylab functionality you need. More information in ipython may be found [here](http://ipython.scipy.org/moin/). 


### Up Next: ###

1.) **Faster calculations:** Currently each diversity measure is calculated from the raw genotypes each time its calculator function is called. This is slow if you wish to calculate more than one diversity measure at a time. Allowing function calls to use pre-calculated data should speed thing up considerably. For example, calculating *Gst* requires calculating *Hs* and *Ht*.  The way the code is written now a subsequent calculation *D* involves re-calculating *Hs* and *Ht*. In the future I'll add code to store the relevant intermediate steps.

2.) **Missing Functionality:** bootstrapping analysis and the Arlequin format parser didn't make this version.

3.) **New Functionality:** One of the main goals of this rewrite is to provide command-line tools.  I'll be producing a command-line version of shortly. I'm also planning on including the *D-est-chao* estimator (eq 13 from Jost, 2008).

### Citations ###

*   Crawford, N. G. (2009). SMOGD: software for the measurement of genetic diversity. Molecular Ecology Resources, 10, 556-557. doi:10.1111/j.1755-0998.2009.02801.x
*   Hedrick, P. W. (2005). A standardized genetic differentiation measure. Evolution; international journal of organic evolution, 59(8), 1633-1638.
*   Jost, L. (2008). Gst and its relatives do not measure differentiation. Molecular Ecology, 17(18), 4015-4026. doi:10.1111/mec.2008.17.issue-18
*   Masatoshi Nei. (1973). Analysis of Gene Diversity in Subdivided Populations. Proceedings of the National Academy of Sciences of the United States of America, 70(12 Pt 1-2), 3321.
*   Nei, M., & Chesser, R. K. (1983). Estimation of fixation indices and gene diversities. Annals of human genetics, 47(Pt 3), 253-259.


