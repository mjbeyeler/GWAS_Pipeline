# Installation
If you only need the GSM extraction script (Scripts/Extract_GSM.ipynb), you only need to install the FaST-LMM part. If you want to use the entire GWAS pipeline, R as described below will also be required.

### FaST-LMM

1) You need a Unix environment (WSL (Windows Subsystem for Linux works too)).
2) Install an Anaconda 2 for Linux distribution: https://www.anaconda.com/download/#linux
3) This pipeline uses FaST-LMM from Microsoft Genomics (https://github.com/MicrosoftGenomics/FaST-LMM). Once Anaconda python is installed, download the full GWAS_Pipeline project, and in the *FaST-LMM* folder type `sudo python setup.py install`. After this the FaST-LMM is fully functional.

#### Warning

FaST-LMM was updated only recently. To make sure to have the latest pysnptools, just type the following in your Anaconda command prompt:

`pip uninstall pysnptools`  
`pip install pysnptool`



### R (for phenotype extraction)
Install any R distribution. Then type

```
sudo apt install libssl-dev        # openssl compatibility
sudo apt install libxml2-dev       # "
sudo apt install gfortran          # "  
```

This (should) make(s) sure that Unix R runs the phenodype adjustments without errors.


# GWAS Pipeline Info

This pipeline

1. adjusts your raw individual phenotype measurements for inversions and *Wolbachia* infections
2. starting from the original DGRP2 variant files (.bed, .bim, .fam triad) from the official DGRP2 website, filters to only keep variants with MAF>0.05. The number of kept variants will be specific to your phenotype, i.e. depending on which lines where measured in your experiments.
3. performs single SNP GWAS on the adjusted phenotype and filtered variants, and stores the results as a file.

   * This pipeline can do permutations.

The following files have to be added to the repositories manually:

* dgrp2.bed, \*.bim and \*.fam to Data/

### Running the Pipeline

Type `./Bash_GWAS_Pipeline_Full.sh` in an unix environment of your choice.

#### Modifying Pipeline Parameters

All important parameters are located in the headers of the files `PipelinePart1_AdjustingPhenotypes_BuildingAlleleFilteringScript.Rmd` and `PipelinePart3_GWASWithPermutations.ipynb`, and can be customized at will.


# GSM Extraction
Open `Scripts/Extract_GSM.ipynb` in Jupyter. More information is given in the notebook itself.