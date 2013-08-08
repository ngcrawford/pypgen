Tutorial:
=========

Once pygen is installed two scripts, ``vcfWindowedFstats`` and ``vcfSNVfstats``, should be available at the command line.

Running [script name].py will print out a short list of commands and adding the ``--help`` or ``-h`` prints out a more detailed list. 

Basic analysis
++++++++++++++

1. Run your samples through GATK or samtools (or similar SNV caller) that emits calls in the standard `VCF format <http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41>`_. By default pypgen's VCF parser only looks at SNVs where the ``FILTER`` column is set to ``PASS`` so you should filter or recalibrate your VCF appropriately before running pypgen.

2. Once you have a VCF file you'll need to bgzip it. `Tabix <http://samtools.sourceforge.net/tabix.shtml>`_ include bgzip so you make sure you have tabix installed. Tabix and samtools installation is detailed in the :ref:`samtools-tabix` section of this guide. The basic command to run bgzip is:

	::
	
 		bgzip -c  path/to/vcf_file.vcf > path/to/vcf_file.vcf.bgz

 This can exceed 30 minutes if your uncompressed VCF file is very large. 

3. Next you need to index your bgzipped VCF file. The command to do this is:

	::
	
		tabix -p vcf path/to/vcf_file.vcf.bgz
	

 This command will produce a ``path/to/vcf_file.vcf.tbi`` index file. 


4. Now you can run pypgen. In a text editor I recommend composing a test command that looks something like this. 

	::
	
	    vcfSNVfstats \
		 -i path/to/vcf_file.vcf.bgz \
		 -p pop1:sample1,sample2 \
		    pop2:sample3,sample4,sample5 \
		    pop3:sample6,sample7,sample8 \
		 -c 2 \
		 -r Chr:1-10001 | head

 You'll need to replace ``path/to/vcf_file.vcf.bgz`` as you did in the last command. 
 
 You'll also need to associate the sample names with their populations. The sample names should to exactly match the sample IDs in the VCF file. If you've forgotten what they are you can run the following command to print them out. 
	
	::
	 
		gunzip -c pypgen/data/example.vcf.gz | grep "#CHROM"
		
 You will also want to change the regions flag such that it selects a valid region in your VCF file. 
 
 Piping the output into ``head`` prevents flooding your terminal with output. 
 
 If you have an enormous number of samples and get an error like ``Argument list too long`` you can just save the text file as a shell script and run it like:
 
 	::
	
		sh path/to/shell_script.sh
		
 If everything worked you should see a header line followed by ~ 9 lines of data. The amount of output varies depending on the region so it's a good idea to pick a region that you know contains SNVs.

 Replacing ``vcfSNVfstats`` with ``vcfWindowedFstats`` and setting the ``--window`` flag is all that is necessary to run a sliding window analysis

Followup analysis
+++++++++++++++++


