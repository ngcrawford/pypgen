Detailed Installation:
======================

Installing pypgen is very straightforward especially if you are familiar with installing python packages. Just follow the instructions below entering the appropriate commands in terminal.

Samtools and tabix:
+++++++++++++++++++

**In OS X:**

#. Install `Xcode <http://itunes.apple.com/us/app/xcode/id497799835>`_ or `Xcode Command Line Tools <https://developer.apple.com/downloads>`_. The CLI tools take up less space, but are an optional install under Xcode. Details on how to do this may be found in the `homebrew documenation <https://github.com/mxcl/homebrew/wiki/Installation#wiki-fn3>`_.

#. Once xcode is installed, install the `homebrew <http://mxcl.github.com/homebrew/>`_ package installer:

	::

		$ ruby -e "$(curl -fsSkL raw.github.com/mxcl/homebrew/go)"

#. Then install samtools using homebrew:

	::
		
		$ brew tap homebrew/science
		$ brew install samtools
		$ brew install tabix

 While you're at it you might want to use brew to install ``wget``.

	:: 

		$ brew install wget

**From source code** (e.g., on linux):

#. Download the latest version of samtools and tabix.

	::

		# replace the ###version### with the appropriate version number (e.g., 0.2.6)
		
		$ wget http://sourceforge.net/projects/samtools/files/tabix/tabix-###version###.tar.bz2
		$ wget http://sourceforge.net/projects/samtools/files/samtools/0.1.18/samtools-0.1.18.tar.bz2

#. Extract the tar.bz2 files

	:: 

		tar jxf *.tar.bz2

#. Then run ``make`` in each directory

#. You'll need to add these directories to your system profile files (e.g., ``.bashrc`` or ``.bash_profile``)

# You can check that everything is working by opening a fresh shell. The commands ``samtools`` and ``tabix`` should now be available from anywhere in the file system.


Pip or Setuptools:
++++++++++++++++++

Note: Pip is the the replacement for setuptools and is the recommended approach.

#. Install pip:

	::
	
		$ curl -O https://raw.github.com/pypa/pip/master/contrib/get-pip.py
		$ [sudo] python get-pip.py
		$ rm get-pip.py

#. But, if you must, you can use setuptools. You'll need to download the appropriate `python .egg file <http://pypi.python.org/pypi/setuptools#files>`_. Then you can run it as an installation script. 

	::
	
		$ [sudo] sh setuptools-0.6c9-py2.4.egg

Pypgen:
+++++++

Pypgen can be installed from `PyPi <http://pypi.python.org/pypi/pypgen>`_  using `pip <http://pypi.python.org/pypi/pip>`_ or `setuptools <http://pypi.python.org/pypi/setuptools>`_:

	::

		$ [sudo] pip install pypgen

or, 

	::

		$ [sudo] easy_install -U pypgen

However, it's recommended, at least in these early days of pypgen when I'm actively fixing bugs, to install it directly from the github repository:


	::

		$ [sudo] pip install -e git+https://github.com/ngcrawford/pypgen.git#egg=Package


This should complete your install. 	   


UnitTests:
++++++++++

UnitTests are in ``pygen/tests/tests.py``

Pysam:
++++++

**NOTE: Pysam should automatically install when you install pypgen. These instructions are only necessary if you have problems with it.**

Pysam is a bit of a finicky installation. The newest versions, in particular seem to have a lot of problems linking to their compiled cython libraries. With that in mind I recommend installing 0.6.

	::

		$ [sudo] pip install pysam==0.6

or, 

	::

		$ [sudo] easy_install -U pysam==0.6

If that doesn't work you'll want to try installing it from source:

	::

		# replace the ###version### with the appropriate version number (e.g., 0.7.4)
		
		$ wget http://pysam.googlecode.com/files/pysam-###version###.tar.gz
		$ tar xzf pysam-###version###.tar.gz
		
then ``cd`` into the directory and run:

	::
	
		$ [sudo] python setup.py install

There is also a `pysam google group <https://groups.google.com/forum/#!forum/pysam-user-group>`_ that is a good source of information.




