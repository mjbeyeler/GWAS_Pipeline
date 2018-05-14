# GWAS Pipeline

The following files have to be added to the repositories manually:
* dgrp2.bed/bim/fam/vcf to Data
* MitoSeq_AllRuns_dm6_chrM.annot.biallellic.vcf to Data
* MitoSeq_AllRuns_dm6_chrM.annot.biallellic_ConvertedReference.bed/bim/fam to Outputs/Plinkfiles.


Other filtered variants of the original dgrp2.vcf or Plink binary files can always be produced at will, depending on your individual needs. :)


## BEWARE

1) In order for FaST-LMM to work, you need Python 2.x. Python 3.x will not work.
2) FaST-LMM was updated only recently. In order to make sure to have the latest pysnptools, just type the following in your Anaconda 2 command prompt.
	* `pip uninstall pysnptools`  
 `pip install pysnptool`