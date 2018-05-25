# Installation
If you only need the GSM extraction script (Scripts/Extract_GSM.ipynb), you only need to install the FaST-LMM part. If you want to use the entire GWAS pipeline, R as described below will also be required.

### FaST-LMM

1. You need a Unix environment. (WSL (Windows Subsystem for Linux) works too.)
2. Install an Anaconda 2 for Linux distribution: https://www.anaconda.com/download/#linux
3. This pipeline uses FaST-LMM from Microsoft Genomics (https://github.com/MicrosoftGenomics/FaST-LMM). Once Anaconda python is installed, download the full GWAS_Pipeline project, and in its *FaST-LMM* folder type `sudo python setup.py install`. After this, FaST-LMM is fully functional.

**Warning**: FaST-LMM was updated only recently. To make sure to have the latest pysnptools, just type the following in your Anaconda command prompt:  
`pip uninstall pysnptools`  
`pip install pysnptool`

### R (for phenotype adjustment)
Install any R distribution in your Unix environment. Then type

```
sudo apt install libssl-dev        # openssl compatibility
sudo apt install libxml2-dev       # "
sudo apt install gfortran          # "  
```

This (should) make(s) sure that Unix R runs the phenodype adjustments without errors.


# GWAS Pipeline Info

This pipeline

1) adjusts your raw individual phenotype measurements for inversions and *Wolbachia* infections
2) Starting with the variants present in the original DGRP2 variant files (.bed, .bim, .fam triad) from the official DGRP2 website, filters to only keep variants with MAF>WhateverThresholdYouChoose (default is 0.05). The number of kept variants will be specific to your phenotype, i.e. depending on which lines where used in the specific experiments
3) performs single SNP GWAS on the adjusted phenotype and filtered variants, and stores the results as a file in the `Outputs/` directory.

   * This pipeline can do permutations.

## Data to be added manually

* `dgrp2.bed`, `*.bim` and `*.fam`, `freeze2.common.rel.mat`, `wolbachia.csv`, `inversion.csv` (all from the official DGRP2 site, `wolbachia` and `inversion` have to be converted from .xlsx to .csv) to `Raw_Data/`
* Your phenotype files to `Inputs/`. They must be named `*Phenotype_Name*_Phenotype_Full.txt`

## Running the Pipeline

1) Put your phenotype in the `Inputs/` folder and name it `*Pheno_Name*_Phenotype_Full.txt`. It should follow the formatting `line_id \t phenotype value \t sex (m/f)`, and should not have a header.

2) In the unix environment of your choice, in the main folder of the project, type  
`var1=value1 var2=value2 ... ./Bash_GWAS_Pipeline_Full.sh\`  
<br/>  Any of the variables are optional and can be omitted! The variables that you can choose from, are:  

Variable | Description
-- | --
pheno | Name of your phenotype (Must correspond ot the way you named it in the phenotype file)
sex | Female / Male / Dimorphism
maf | minor allele frequency, between 0 and 0.5 (default is 0.05)
perm | the number of permutations you want to perform (default is 0)
reproducibility | If you want to run the script with all package versions in the exact state as when the scripts were written, put this to TRUE (default is TRUE)
use_official_gsm | If set to FALSE, will calculate GSM from provided variants (default). If set to TRUE, will use the Freeze 2 GSM provided by Mackay group.
  
**Example:** `pheno=Mass sex=Male perm=5 reproducibility=FALSE ./Bash_GWAS_Pipeline_Full.sh`

# GSM Extraction
Open `Scripts/Extract_GSM.ipynb` in Jupyter. More information is given in the notebook itself.