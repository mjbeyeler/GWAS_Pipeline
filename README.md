# GWAS Pipeline Info

This pipeline

1. adjusts your raw individual phenotype measurements for inversions and *Wolbachia* infections
2. starting from the original DGRP2 variant files (.bed, .bim, .fam triad) from the official DGRP2 website, filters to only keep variants with MAF>0.05. The number of kept variants will be specific to your phenotype, i.e. depending on which lines where measured in your experiments.
3. performs single SNP GWAS on the adjusted phenotype and filtered variants, and stores the results as a file.

   * This pipeline can do permutations.

The following files have to be added to the repositories manually:

* dgrp2.bed, \*.bim and \*.fam to Data/

# Running the Pipeline

Type `./Bash_GWAS_Pipeline_Full.sh` in an unix environment of your choice.

First run these commands to make sure Unix R will work without any errors:
```
sudo apt install libssl-dev        # openssl  
sudo apt install libxml2-dev       # to be able to install igraph package  
sudo apt install gfortran          # igraph compatibility  
```

# FaST-LMM Warning

1. In order for FaST-LMM to work, you need Python 2.x. Python 3.x will not work.
2. FaST-LMM was updated only recently. In order to make sure to have the latest pysnptools, just type the following in your Anaconda 2 command prompt.
* `pip uninstall pysnptools`  
  `pip install pysnptool`