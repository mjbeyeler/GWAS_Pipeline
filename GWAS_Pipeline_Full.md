GWAS Pipeline Full
================
Michael Beyeler
2018-05-13

Setup
=====

Loading packages, Reproducibility
---------------------------------

**Concerning reproducibility**: In order to guarantee reproducibility, keep the `.BuildReproducibleEnvironment(...)` FUNCTION active.

    ## Warning: Your R version (3.5.0) differs from the version this project was programmed in (3.4.3),
    ## which might be causing invonveniences.

    ## Warning in library(package, lib.loc = lib.loc, character.only = TRUE,
    ## logical.return = TRUE, : there is no package called 'data.table'

    ## package 'data.table' successfully unpacked and MD5 sums checked
    ## 
    ## The downloaded binary packages are in
    ##  C:\Users\micha\AppData\Local\Temp\RtmpmWditr\downloaded_packages

    ## Warning in library(package, lib.loc = lib.loc, character.only = TRUE,
    ## logical.return = TRUE, : there is no package called 'icesTAF'

    ## package 'jsonlite' successfully unpacked and MD5 sums checked
    ## package 'mime' successfully unpacked and MD5 sums checked
    ## package 'curl' successfully unpacked and MD5 sums checked
    ## package 'openssl' successfully unpacked and MD5 sums checked
    ## package 'R6' successfully unpacked and MD5 sums checked
    ## package 'httr' successfully unpacked and MD5 sums checked
    ## package 'icesTAF' successfully unpacked and MD5 sums checked
    ## 
    ## The downloaded binary packages are in
    ##  C:\Users\micha\AppData\Local\Temp\RtmpmWditr\downloaded_packages

    ## Warning in library(package, lib.loc = lib.loc, character.only = TRUE,
    ## logical.return = TRUE, : there is no package called 'lintr'

    ## 
    ##   There is a binary version available but the source version is
    ##   later:
    ##         binary source needs_compilation
    ## stringi  1.1.7  1.2.2              TRUE
    ## 
    ##   Binaries will be installed
    ## package 'assertthat' successfully unpacked and MD5 sums checked
    ## package 'glue' successfully unpacked and MD5 sums checked
    ## package 'stringi' successfully unpacked and MD5 sums checked
    ## package 'magrittr' successfully unpacked and MD5 sums checked
    ## package 'lazyeval' successfully unpacked and MD5 sums checked
    ## package 'cli' successfully unpacked and MD5 sums checked
    ## package 'praise' successfully unpacked and MD5 sums checked
    ## package 'rlang' successfully unpacked and MD5 sums checked
    ## package 'withr' successfully unpacked and MD5 sums checked
    ## package 'pkgconfig' successfully unpacked and MD5 sums checked
    ## package 'evaluate' successfully unpacked and MD5 sums checked
    ## package 'highr' successfully unpacked and MD5 sums checked
    ## package 'markdown' successfully unpacked and MD5 sums checked
    ## package 'stringr' successfully unpacked and MD5 sums checked
    ## package 'yaml' successfully unpacked and MD5 sums checked
    ## package 'rex' successfully unpacked and MD5 sums checked
    ## package 'crayon' successfully unpacked and MD5 sums checked
    ## package 'stringdist' successfully unpacked and MD5 sums checked
    ## package 'testthat' successfully unpacked and MD5 sums checked
    ## package 'digest' successfully unpacked and MD5 sums checked
    ## package 'igraph' successfully unpacked and MD5 sums checked
    ## package 'rstudioapi' successfully unpacked and MD5 sums checked
    ## package 'knitr' successfully unpacked and MD5 sums checked
    ## package 'lintr' successfully unpacked and MD5 sums checked
    ## 
    ## The downloaded binary packages are in
    ##  C:\Users\micha\AppData\Local\Temp\RtmpmWditr\downloaded_packages

    ## Warning in library(package, lib.loc = lib.loc, character.only = TRUE,
    ## logical.return = TRUE, : there is no package called 'lme4'

    ## package 'minqa' successfully unpacked and MD5 sums checked
    ## package 'nloptr' successfully unpacked and MD5 sums checked
    ## package 'Rcpp' successfully unpacked and MD5 sums checked
    ## package 'RcppEigen' successfully unpacked and MD5 sums checked
    ## package 'lme4' successfully unpacked and MD5 sums checked
    ## 
    ## The downloaded binary packages are in
    ##  C:\Users\micha\AppData\Local\Temp\RtmpmWditr\downloaded_packages

    ## Warning in library(package, lib.loc = lib.loc, character.only = TRUE,
    ## logical.return = TRUE, : there is no package called 'reticulate'

    ## package 'reticulate' successfully unpacked and MD5 sums checked
    ## 
    ## The downloaded binary packages are in
    ##  C:\Users\micha\AppData\Local\Temp\RtmpmWditr\downloaded_packages

    ## Warning in library(package, lib.loc = lib.loc, character.only = TRUE,
    ## logical.return = TRUE, : there is no package called 'tictoc'

    ## package 'tictoc' successfully unpacked and MD5 sums checked
    ## 
    ## The downloaded binary packages are in
    ##  C:\Users\micha\AppData\Local\Temp\RtmpmWditr\downloaded_packages

``` r
cat("The current checkpoint is:\n")
```

    ## The current checkpoint is:

``` r
checkpoint::setSnapshot()
```

    ##                                             CRAN 
    ## "https://mran.microsoft.com/snapshot/2018-05-15"

``` r
cat("\n\nSession Info:\n\n")
```

    ## 
    ## 
    ## Session Info:

``` r
sessionInfo()
```

    ## R version 3.5.0 (2018-04-23)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 17134)
    ## 
    ## Matrix products: default
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_United States.1252 
    ## [2] LC_CTYPE=English_United States.1252   
    ## [3] LC_MONETARY=English_United States.1252
    ## [4] LC_NUMERIC=C                          
    ## [5] LC_TIME=English_United States.1252    
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] tictoc_1.0        reticulate_1.7    lme4_1.1-17       Matrix_1.2-14    
    ## [5] lintr_1.0.2       icesTAF_1.5-3     data.table_1.11.2 checkpoint_0.4.3 
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] rex_1.1.2       Rcpp_0.12.16    knitr_1.20      magrittr_1.5   
    ##  [5] MASS_7.3-49     splines_3.5.0   lattice_0.20-35 R6_2.2.2       
    ##  [9] minqa_1.2.4     stringr_1.3.1   httr_1.3.1      tools_3.5.0    
    ## [13] grid_3.5.0      nlme_3.1-137    htmltools_0.3.6 yaml_2.1.19    
    ## [17] lazyeval_0.2.1  rprojroot_1.3-2 digest_0.6.15   nloptr_1.0.4   
    ## [21] evaluate_0.10.1 rmarkdown_1.9   stringi_1.1.7   compiler_3.5.0 
    ## [25] backports_1.1.2 jsonlite_1.5

Constants
=========

``` r
SEXUAL.DIMORPHISM <- TRUE
PHENOTYPE.NAME <- 'Food_Intake'
NORMALITY.SIGNIFICANCE.LEVEL <- 0.05
INVERSIONS.CONSIDERED <- c('In.2L.t', 'In.2R.NS', 'In.3R.P', 'In.3R.K', 'In.3R.Mo')

MAF.THRESHOLD=0.25
```

``` python
NUMBER_OF_PERMUTATIONS = 0
OUTPUT_NAME = '../Outputs/GWASMay'
phenotype_data = '../Outputs/Fast-Lmm-Input-Mass-Female.txt'
variants_to_test = '../Outputs/Current_Pipeline_Variants'
```

Required Files
==============

Data
----

``` r
# Main data

Phenotype_Raw <- read.delim('Inputs/Mass-Data-with-Line-IDs.txt', header=T, sep=" ")
# Phenotype_Raw <- read.delim('Data/Food-Intake-Garlapow.csv', header=F, sep=',')
# Phenotype_Raw <- read.delim('Data/Vonesch2016-IOD-Raw.txt', header=F)

# Supporting data

Dgrp2_Inversions <- read.csv('Data/inversion.csv', header=T)
Dgrp2_Infection <- read.csv('Data/wolbachia.csv')

# Some phenotypes actually use the flystock ID instead of the DGRP ID.
# For these cases, this data frame will come in handy.
# Dgrp_Flystock_Ids <- read.delim('Data/Dgrp-Flystocks-Ids.txt',
#                                 comment.char='#')
```

Functions
---------

Functions that were specifically programmed for this script

``` r
source('Functions/NormalityHistogram.R')
source('Functions/ChisqForNormality.R')
source('Functions/WriteBare.R')
```

Setting up Anaconda python
--------------------------

Making sure that Anaconda 2 python is used in case multiple python distributions are installed. FaST-LMM only works on Anaconda 2 python.

Path is set up depending on what OS is used

To-do: check if it also works for Mac

``` r
if(.Platform$OS.type == "unix") {
  knitr::opts_chunk$set(engine.path = list(python = '/anaconda/bin/python'))
  # use_python('/anaconda/bin/python')
} else {
  knitr::opts_chunk$set(engine.path = list(python = file.path(Sys.getenv("USERPROFILE"),"Anaconda2\\python.exe", fsep='\\')))
  # use_python(file.path(Sys.getenv("USERPROFILE"),"Anaconda2\\python.exe", fsep="\\"))
}
```

Adjustment procedure
--------------------

    ## [1] "Unique to Female:  "
    ## [1] "Unique to Male:  "

![](GWAS_Pipeline_Full_files/figure-markdown_github/unnamed-chunk-5-1.png)![](GWAS_Pipeline_Full_files/figure-markdown_github/unnamed-chunk-5-2.png)![](GWAS_Pipeline_Full_files/figure-markdown_github/unnamed-chunk-5-3.png)![](GWAS_Pipeline_Full_files/figure-markdown_github/unnamed-chunk-5-4.png)![](GWAS_Pipeline_Full_files/figure-markdown_github/unnamed-chunk-5-5.png)

    ## [1] "p-value total 5.2630561996935e-21"
    ## [1] 1.025024e-05
    ## [1] 0.008176191
    ## [1] "At least one of the phenotypes was not normally distributed. Log-transformation performed."

![](GWAS_Pipeline_Full_files/figure-markdown_github/unnamed-chunk-5-6.png)![](GWAS_Pipeline_Full_files/figure-markdown_github/unnamed-chunk-5-7.png)

Filtering Lines and Low Minor Allele Frequencies
================================================

At the moment, bash is causing some trouble in Windows R Markdown. Thus, in the meantime, I'm using a workaround to write a bash shell script using R code, and then executing it in an unix environment.

``` r
cat(paste("#!/bin/bash

PHENOTYPE_NAME=", PHENOTYPE.NAME,"
MAF=", MAF.THRESHOLD,"

cd plink2_linux_x86_64
./plink2 --bfile ../Data/dgrp2 --keep ../Outputs/Plink-Lines-$PHENOTYPE_NAME.txt --maf $MAF --make-bed --out ../Outputs/Current_Pipeline_Variants
", sep=''),
file='Scripts/Plink2_Filtering_Alleles.sh')

# The following command is necessary to 
dos2unix('Scripts/Plink2_Filtering_Alleles.sh')
```

``` bash
./Scripts/Plink2_Filtering_Alleles.sh
```

    ## PLINK v2.00a1LM 64-bit Intel (11 Feb 2018)     www.cog-genomics.org/plink/2.0/
    ## (C) 2005-2018 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to ../Outputs/Current_Pipeline_Variants.log.
    ## Options in effect:
    ##   --bfile ../Data/dgrp2
    ##   --keep ../Outputs/Plink-Lines-Food_Intake.txt
    ##   --maf 0.25
    ##   --make-bed
    ##   --out ../Outputs/Current_Pipeline_Variants
    ## 
    ## Start time: Tue May 15 20:45:59 2018
    ## 16221 MB RAM detected; reserving 8110 MB for main workspace.
    ## Using up to 8 compute threads.
    ## 205 samples (205 females, 0 males; 205 founders) loaded from ../Data/dgrp2.fam.
    ## 4438427 variants loaded from ../Data/dgrp2.bim.
    ## Note: No phenotype data present.
    ## --keep: 157 samples remaining.
    ## 157 samples (157 females, 0 males; 157 founders) remaining after main filters.
    ## Calculating allele frequencies... 0%1%2%4%5%7%8%10%11%13%14%16%17%19%20%22%23%25%26%28%29%31%32%33%35%36%38%39%41%42%44%45%47%48%50%51%53%54%56%57%59%60%62%63%64%66%67%69%70%72%73%75%76%78%79%81%82%84%85%87%88%90%91%93%94%95%97%98%done.
    ## 3682349 variants removed due to minor allele threshold(s)
    ## (--maf/--max-maf/--mac/--max-mac).
    ## 756078 variants remaining after main filters.
    ## Writing ../Outputs/Current_Pipeline_Variants.bed ... 0%1%3%4%5%7%9%10%12%13%15%16%18%19%20%22%23%24%26%27%29%31%32%34%35%37%38%40%41%43%44%46%47%49%50%52%54%55%57%58%60%61%63%64%65%67%68%69%71%72%74%76%77%78%80%81%83%84%86%87%89%90%92%93%95%96%98%done.
    ## Writing ../Outputs/Current_Pipeline_Variants.bim ... done.
    ## Writing ../Outputs/Current_Pipeline_Variants.fam ... done.
    ## End time: Tue May 15 20:46:00 2018

FaST-LMM GWAS
=============

``` python
# coding: utf-8
# # EXPLANATION
# This script performs GWAS on a pre-selected number of variants (in Plink binary format) and on a specific phenotype you want. Ideally, the phenotype has been adjusted for inversions, infection, and also log-transformed if the raw phenotype wasn't deemed normally distributed.
#
# BEWARE: You should format your phenotype line IDs in the same way as denoted in the .fam variant file you are using.
#  # Time measurement
#
# In[1]:
import time
start = time.clock()
# # Importing, general preparations
# In[2]:
import os
os.chdir('FaST_LMM')
from fastlmm.association import single_snp
```

    ## C:\Users\micha\ANACON~1\lib\site-packages\h5py\__init__.py:36: FutureWarning: Conversion of the second argument of issubdtype from `float` to `np.floating` is deprecated. In future, it will be treated as `np.float64 == np.dtype(float).type`.
    ##   from ._conv import register_converters as _register_converters
    ## C:\Users\micha\ANACON~1\lib\site-packages\sklearn\cross_validation.py:41: DeprecationWarning: This module was deprecated in version 0.18 in favor of the model_selection module into which all the refactored classes and functions are moved. Also note that the interface of the new CV iterators are different from that of this module. This module will be removed in 0.20.
    ##   "This module will be removed in 0.20.", DeprecationWarning)

``` python
from fastlmm.inference.fastlmm_predictor import _snps_fixup, _pheno_fixup, _kernel_fixup, _SnpTrainTest
from random import shuffle
import numpy as np
import pandas as pd
import time
import os
# We're going to need PySnpTools, to do permutations, because we can shuffle bed files by varint using Bed
import sys
sys.path.append('../PySnpTools')
from pysnptools.snpreader import Bed
from shutil import copyfile
# # CONSTANTS
# In[3]:
# Permutation info
#NUMBER_OF_PERMUTATIONS = 0
# Where you want to save your phenotype
#OUTPUT_NAME = '../Outputs/GWASMay'
# # INPUTS
# In[4]:
# Phenotype data to test
#phenotype_data = '../Outputs/Fast-Lmm-Input-Interocular-Distance-Vonesch2016-Female.txt'
#phenotype_data = '../Outputs/Fast-Lmm-Input-Mass-Female.txt'
#phenotype_data = '../Outputs/Dgrp2-AllLines-RandomPheno-for-Mito.txt'
#phenotype_data = '../Outputs/T5-Pheno-for-Fast-Lmm-DGRPFormat.txt'
#phenotype_data = '../Outputs/T5-Pheno-for-Fast-Lmm.txt'
# Variants to test your phenotype on
#variants_to_test = '../Outputs/MitoSeq_AllRuns_dm6_chrM.annot.biallellic_ConvertedReference'
#variants_to_test = '../Data/dgrp2'
#variants_to_test = '../Outputs/Plinkfiles/DGRP2Chr2R'
#variants_to_test = '../Outputs/Current_Pipeline_Variants'
#variants_to_test = '../Outputs/Plinkfiles/Dgrp2-CSLines-Maf5'
# ACTUAL GWAS =================================================================
# Clearing cache
# This ensures that the relationship matrix is recalculated for each phenotype.
try:
    os.remove('Outputs/Fast-Lmm-Cache/Gwas-Permutations-Cache.npz')
except OSError:
    pass
# Performing GWAS on the real phenotype:
time_0 = time.time()
results_df = single_snp(variants_to_test,  phenotype_data,
                        cache_file='Outputs/Fast-Lmm-Cache/Gwas-Permutations-Cache.npz',
                        leave_out_one_chrom=False,
                        save_test_statistic=True,
                        output_file_name = OUTPUT_NAME + '-Original.txt',
                        )
```

    ## C:\Users\micha\ANACON~1\lib\site-packages\pysnptools-0.3.13-py2.7-win-amd64.egg\pysnptools\snpreader\bed.py:42: FutureWarning: 'count_A1' was not set. For now it will default to 'False', but in the future it will default to 'True'
    ##   warnings.warn("'count_A1' was not set. For now it will default to 'False', but in the future it will default to 'True'", FutureWarning)
    ## C:\Users\micha\ANACON~1\lib\site-packages\pysnptools-0.3.13-py2.7-win-amd64.egg\pysnptools\snpreader\snpreader.py:625: FutureWarning: Conversion of the second argument of issubdtype from `str` to `str` is deprecated. In future, it will be treated as `np.string_ == np.dtype(str).type`.
    ##   assert np.issubdtype(self._row.dtype, str) and len(self._row.shape)==2 and self._row.shape[1]==2, "iid should be dtype str, have two dimensions, and the second dimension should be size 2"
    ## C:\Users\micha\ANACON~1\lib\site-packages\pysnptools-0.3.13-py2.7-win-amd64.egg\pysnptools\snpreader\snpreader.py:626: FutureWarning: Conversion of the second argument of issubdtype from `str` to `str` is deprecated. In future, it will be treated as `np.string_ == np.dtype(str).type`.
    ##   assert np.issubdtype(self._col.dtype, str) and len(self._col.shape)==1, "sid should be of dtype of str and one dimensional"
    ## C:\Users\micha\ANACON~1\lib\site-packages\pysnptools-0.3.13-py2.7-win-amd64.egg\pysnptools\kernelreader\kernelreader.py:325: FutureWarning: Conversion of the second argument of issubdtype from `str` to `str` is deprecated. In future, it will be treated as `np.string_ == np.dtype(str).type`.
    ##   assert np.issubdtype(self._row.dtype, str) and len(self._row.shape)==2 and self._row.shape[1]==2, "iid0 should be dtype str, have two dimensions, and the second dimension should be size 2"
    ## C:\Users\micha\ANACON~1\lib\site-packages\pysnptools-0.3.13-py2.7-win-amd64.egg\pysnptools\kernelreader\kernelreader.py:326: FutureWarning: Conversion of the second argument of issubdtype from `str` to `str` is deprecated. In future, it will be treated as `np.string_ == np.dtype(str).type`.
    ##   assert np.issubdtype(self._col.dtype, str) and len(self._col.shape)==2 and self._col.shape[1]==2, "iid1 should be dtype str, have two dimensions, and the second dimension should be size 2"
    ## C:\Users\micha\ANACON~1\lib\site-packages\fastlmm-0.2.31-py2.7-win-amd64.egg\fastlmm\inference\lmm_cov.py:793: RuntimeWarning: invalid value encountered in divide
    ##   beta = snpsKY / (snpsKsnps + penalty_)
    ## C:\Users\micha\ANACON~1\lib\site-packages\numpy\core\_methods.py:29: RuntimeWarning: invalid value encountered in reduce
    ##   return umr_minimum(a, axis, None, out, keepdims)
    ## Traceback (most recent call last):
    ##   File "C:\Users\micha\ANACON~1\lib\logging\__init__.py", line 882, in emit
    ##     stream.write(fs % msg)
    ##   File "C:\Users\micha\ANACON~1\lib\StringIO.py", line 213, in write
    ##     _complain_ifclosed(self.closed)
    ##   File "C:\Users\micha\ANACON~1\lib\StringIO.py", line 40, in _complain_ifclosed
    ##     raise ValueError, "I/O operation on closed file"
    ## ValueError: I/O operation on closed file
    ## Logged from file lmm_cov.py, line 795
    ## C:\Users\micha\ANACON~1\lib\site-packages\fastlmm-0.2.31-py2.7-win-amd64.egg\fastlmm\inference\lmm_cov.py:804: RuntimeWarning: divide by zero encountered in divide
    ##   variance_beta = r2 / (N - 1.0) / snpsKsnps
    ## Traceback (most recent call last):
    ##   File "C:\Users\micha\ANACON~1\lib\logging\__init__.py", line 882, in emit
    ##     stream.write(fs % msg)
    ##   File "C:\Users\micha\ANACON~1\lib\StringIO.py", line 213, in write
    ##     _complain_ifclosed(self.closed)
    ##   File "C:\Users\micha\ANACON~1\lib\StringIO.py", line 40, in _complain_ifclosed
    ##     raise ValueError, "I/O operation on closed file"
    ## ValueError: I/O operation on closed file
    ## Logged from file lmm_cov.py, line 795
    ## Traceback (most recent call last):
    ##   File "C:\Users\micha\ANACON~1\lib\logging\__init__.py", line 882, in emit
    ##     stream.write(fs % msg)
    ##   File "C:\Users\micha\ANACON~1\lib\StringIO.py", line 213, in write
    ##     _complain_ifclosed(self.closed)
    ##   File "C:\Users\micha\ANACON~1\lib\StringIO.py", line 40, in _complain_ifclosed
    ##     raise ValueError, "I/O operation on closed file"
    ## ValueError: I/O operation on closed file
    ## Logged from file lmm_cov.py, line 795
    ## Traceback (most recent call last):
    ##   File "C:\Users\micha\ANACON~1\lib\logging\__init__.py", line 882, in emit
    ##     stream.write(fs % msg)
    ##   File "C:\Users\micha\ANACON~1\lib\StringIO.py", line 213, in write
    ##     _complain_ifclosed(self.closed)
    ##   File "C:\Users\micha\ANACON~1\lib\StringIO.py", line 40, in _complain_ifclosed
    ##     raise ValueError, "I/O operation on closed file"
    ## ValueError: I/O operation on closed file
    ## Logged from file lmm_cov.py, line 795
    ## Traceback (most recent call last):
    ##   File "C:\Users\micha\ANACON~1\lib\logging\__init__.py", line 882, in emit
    ##     stream.write(fs % msg)
    ##   File "C:\Users\micha\ANACON~1\lib\StringIO.py", line 213, in write
    ##     _complain_ifclosed(self.closed)
    ##   File "C:\Users\micha\ANACON~1\lib\StringIO.py", line 40, in _complain_ifclosed
    ##     raise ValueError, "I/O operation on closed file"
    ## ValueError: I/O operation on closed file
    ## Logged from file lmm_cov.py, line 795
    ## Traceback (most recent call last):
    ##   File "C:\Users\micha\ANACON~1\lib\logging\__init__.py", line 882, in emit
    ##     stream.write(fs % msg)
    ##   File "C:\Users\micha\ANACON~1\lib\StringIO.py", line 213, in write
    ##     _complain_ifclosed(self.closed)
    ##   File "C:\Users\micha\ANACON~1\lib\StringIO.py", line 40, in _complain_ifclosed
    ##     raise ValueError, "I/O operation on closed file"
    ## ValueError: I/O operation on closed file
    ## Logged from file lmm_cov.py, line 795
    ## Traceback (most recent call last):
    ##   File "C:\Users\micha\ANACON~1\lib\logging\__init__.py", line 882, in emit
    ##     stream.write(fs % msg)
    ##   File "C:\Users\micha\ANACON~1\lib\StringIO.py", line 213, in write
    ##     _complain_ifclosed(self.closed)
    ##   File "C:\Users\micha\ANACON~1\lib\StringIO.py", line 40, in _complain_ifclosed
    ##     raise ValueError, "I/O operation on closed file"
    ## ValueError: I/O operation on closed file
    ## Logged from file lmm_cov.py, line 795
    ## Traceback (most recent call last):
    ##   File "C:\Users\micha\ANACON~1\lib\logging\__init__.py", line 882, in emit
    ##     stream.write(fs % msg)
    ##   File "C:\Users\micha\ANACON~1\lib\StringIO.py", line 213, in write
    ##     _complain_ifclosed(self.closed)
    ##   File "C:\Users\micha\ANACON~1\lib\StringIO.py", line 40, in _complain_ifclosed
    ##     raise ValueError, "I/O operation on closed file"
    ## ValueError: I/O operation on closed file
    ## Logged from file lmm_cov.py, line 795
    ## Traceback (most recent call last):
    ##   File "C:\Users\micha\ANACON~1\lib\logging\__init__.py", line 882, in emit
    ##     stream.write(fs % msg)
    ##   File "C:\Users\micha\ANACON~1\lib\StringIO.py", line 213, in write
    ##     _complain_ifclosed(self.closed)
    ##   File "C:\Users\micha\ANACON~1\lib\StringIO.py", line 40, in _complain_ifclosed
    ##     raise ValueError, "I/O operation on closed file"
    ## ValueError: I/O operation on closed file
    ## Logged from file lmm_cov.py, line 795
    ## Traceback (most recent call last):
    ##   File "C:\Users\micha\ANACON~1\lib\logging\__init__.py", line 882, in emit
    ##     stream.write(fs % msg)
    ##   File "C:\Users\micha\ANACON~1\lib\StringIO.py", line 213, in write
    ##     _complain_ifclosed(self.closed)
    ##   File "C:\Users\micha\ANACON~1\lib\StringIO.py", line 40, in _complain_ifclosed
    ##     raise ValueError, "I/O operation on closed file"
    ## ValueError: I/O operation on closed file
    ## Logged from file lmm_cov.py, line 795
    ## Traceback (most recent call last):
    ##   File "C:\Users\micha\ANACON~1\lib\logging\__init__.py", line 882, in emit
    ##     stream.write(fs % msg)
    ##   File "C:\Users\micha\ANACON~1\lib\StringIO.py", line 213, in write
    ##     _complain_ifclosed(self.closed)
    ##   File "C:\Users\micha\ANACON~1\lib\StringIO.py", line 40, in _complain_ifclosed
    ##     raise ValueError, "I/O operation on closed file"
    ## ValueError: I/O operation on closed file
    ## Logged from file lmm_cov.py, line 795
    ## Traceback (most recent call last):
    ##   File "C:\Users\micha\ANACON~1\lib\logging\__init__.py", line 882, in emit
    ##     stream.write(fs % msg)
    ##   File "C:\Users\micha\ANACON~1\lib\StringIO.py", line 213, in write
    ##     _complain_ifclosed(self.closed)
    ##   File "C:\Users\micha\ANACON~1\lib\StringIO.py", line 40, in _complain_ifclosed
    ##     raise ValueError, "I/O operation on closed file"
    ## ValueError: I/O operation on closed file
    ## Logged from file lmm_cov.py, line 795
    ## Traceback (most recent call last):
    ##   File "C:\Users\micha\ANACON~1\lib\logging\__init__.py", line 882, in emit
    ##     stream.write(fs % msg)
    ##   File "C:\Users\micha\ANACON~1\lib\StringIO.py", line 213, in write
    ##     _complain_ifclosed(self.closed)
    ##   File "C:\Users\micha\ANACON~1\lib\StringIO.py", line 40, in _complain_ifclosed
    ##     raise ValueError, "I/O operation on closed file"
    ## ValueError: I/O operation on closed file
    ## Logged from file lmm_cov.py, line 795
    ## Traceback (most recent call last):
    ##   File "C:\Users\micha\ANACON~1\lib\logging\__init__.py", line 882, in emit
    ##     stream.write(fs % msg)
    ##   File "C:\Users\micha\ANACON~1\lib\StringIO.py", line 213, in write
    ##     _complain_ifclosed(self.closed)
    ##   File "C:\Users\micha\ANACON~1\lib\StringIO.py", line 40, in _complain_ifclosed
    ##     raise ValueError, "I/O operation on closed file"
    ## ValueError: I/O operation on closed file
    ## Logged from file lmm_cov.py, line 795
    ## Traceback (most recent call last):
    ##   File "C:\Users\micha\ANACON~1\lib\logging\__init__.py", line 882, in emit
    ##     stream.write(fs % msg)
    ##   File "C:\Users\micha\ANACON~1\lib\StringIO.py", line 213, in write
    ##     _complain_ifclosed(self.closed)
    ##   File "C:\Users\micha\ANACON~1\lib\StringIO.py", line 40, in _complain_ifclosed
    ##     raise ValueError, "I/O operation on closed file"
    ## ValueError: I/O operation on closed file
    ## Logged from file lmm_cov.py, line 795
    ## Traceback (most recent call last):
    ##   File "C:\Users\micha\ANACON~1\lib\logging\__init__.py", line 882, in emit
    ##     stream.write(fs % msg)
    ##   File "C:\Users\micha\ANACON~1\lib\StringIO.py", line 213, in write
    ##     _complain_ifclosed(self.closed)
    ##   File "C:\Users\micha\ANACON~1\lib\StringIO.py", line 40, in _complain_ifclosed
    ##     raise ValueError, "I/O operation on closed file"
    ## ValueError: I/O operation on closed file
    ## Logged from file lmm_cov.py, line 795
    ## Traceback (most recent call last):
    ##   File "C:\Users\micha\ANACON~1\lib\logging\__init__.py", line 882, in emit
    ##     stream.write(fs % msg)
    ##   File "C:\Users\micha\ANACON~1\lib\StringIO.py", line 213, in write
    ##     _complain_ifclosed(self.closed)
    ##   File "C:\Users\micha\ANACON~1\lib\StringIO.py", line 40, in _complain_ifclosed
    ##     raise ValueError, "I/O operation on closed file"
    ## ValueError: I/O operation on closed file
    ## Logged from file lmm_cov.py, line 795
    ## Traceback (most recent call last):
    ##   File "C:\Users\micha\ANACON~1\lib\logging\__init__.py", line 882, in emit
    ##     stream.write(fs % msg)
    ##   File "C:\Users\micha\ANACON~1\lib\StringIO.py", line 213, in write
    ##     _complain_ifclosed(self.closed)
    ##   File "C:\Users\micha\ANACON~1\lib\StringIO.py", line 40, in _complain_ifclosed
    ##     raise ValueError, "I/O operation on closed file"
    ## ValueError: I/O operation on closed file
    ## Logged from file lmm_cov.py, line 795
    ## Traceback (most recent call last):
    ##   File "C:\Users\micha\ANACON~1\lib\logging\__init__.py", line 882, in emit
    ##     stream.write(fs % msg)
    ##   File "C:\Users\micha\ANACON~1\lib\StringIO.py", line 213, in write
    ##     _complain_ifclosed(self.closed)
    ##   File "C:\Users\micha\ANACON~1\lib\StringIO.py", line 40, in _complain_ifclosed
    ##     raise ValueError, "I/O operation on closed file"
    ## ValueError: I/O operation on closed file
    ## Logged from file lmm_cov.py, line 795
    ## Traceback (most recent call last):
    ##   File "C:\Users\micha\ANACON~1\lib\logging\__init__.py", line 882, in emit
    ##     stream.write(fs % msg)
    ##   File "C:\Users\micha\ANACON~1\lib\StringIO.py", line 213, in write
    ##     _complain_ifclosed(self.closed)
    ##   File "C:\Users\micha\ANACON~1\lib\StringIO.py", line 40, in _complain_ifclosed
    ##     raise ValueError, "I/O operation on closed file"
    ## ValueError: I/O operation on closed file
    ## Logged from file lmm_cov.py, line 795

``` python
time_1 = time.time()
print('Time for full GWAS:' + str(time_1 - time_0) + 's')
```

    ## Time for full GWAS:167.193000078s

``` python
print 'Total time: ' + str(time.clock()-start)
# In[ ]:
```

    ## Total time: 167.737829691

``` python
test_stat = pd.read_csv('../Outputs/Fast-Lmm-Cache/Test-Stat-Cache.txt', header=None)
test_stat = test_stat.replace('[\[\] ]', '', regex=True)
test_stat = pd.to_numeric(test_stat[0])
results_df['Full ID'] = results_df['Chr'].astype('str') + '_' + results_df['ChrPos'].astype('str')
results_df = pd.concat([results_df[['Chr', 'ChrPos', 'SNP', 'Full ID', 'PValue']], test_stat],
                       axis = 1)
results_df.columns = ['Chr', 'ChrPos', 'SNP', 'Full ID', 'PValue', 'F-test statistic']
mybed = Bed(variants_to_test + '.bed')
mysnpdata = mybed.read()
print 'Time: ' + str(time.clock()-start)
# In[ ]:
```

    ## Time: 178.764795037

``` python
pheno = _pheno_fixup(phenotype_data, count_A1=None).read()
pheno = pheno.val[np.searchsorted(pheno.iid[:,1], mysnpdata.iid[:,1])]
snpdata = mysnpdata.val
diff = range(snpdata.shape[1])
maf = range(snpdata.shape[1])
n_alleles = range(snpdata.shape[1])
mean_major = range(snpdata.shape[1])
for i in range(snpdata.shape[1]):
    ref = [j for j, x in enumerate(snpdata[:,i]) if x == 2]
    alt = [j for j, x in enumerate(snpdata[:,i]) if x == 0]
    meanref = np.mean(pheno[ref])
    meanalt = np.mean(pheno[alt])
    if len(ref) > len(alt):
        diff[i] = meanref - meanalt
        maf[i] = float(len(alt)) / (len(ref) + len(alt))
        n_alleles[i] = len(ref) + len(alt)
        mean_major[i] = meanref
    elif len(ref) + len(alt) == 0:
        diff[i] = float('NaN')
        maf[i] = float('NaN')
        n_alleles[i] = len(ref) + len(alt)
        mean_major[i] = float('NaN')
    else:
        diff[i] = meanalt - meanref
        maf[i] = float(len(ref)) / (len(ref) + len(alt))
        n_alleles[i] = len(ref) + len(alt)
        mean_major[i] = meanalt
```

    ## C:\Users\micha\ANACON~1\lib\site-packages\numpy\core\fromnumeric.py:2957: RuntimeWarning: Mean of empty slice.
    ##   out=out, **kwargs)
    ## C:\Users\micha\ANACON~1\lib\site-packages\numpy\core\_methods.py:80: RuntimeWarning: invalid value encountered in double_scalars
    ##   ret = ret.dtype.type(ret / rcount)

``` python
print 'Time: ' + str(time.clock()-start)
# In[ ]:
```

    ## Time: 265.108094113

``` python
diff_df = diff_df = pd.DataFrame(data={'MajMinDiff':diff,
                                       'MeanMajor': mean_major,
                                       'MAF':maf,
                                       'NAlleles':n_alleles})
diff_df['SNP'] = mysnpdata.sid
results_df = pd.merge(results_df, diff_df, on='SNP')
```

PHENOTYPE shuffling/permutation and adding the p-values from the resulting
==========================================================================

GWAS to the results data frame.
===============================

phenotype\_to\_shuffle = pd.read\_table(phenotype\_data,
========================================================

sep=' ', header=None)
=====================

indices = range(len(phenotype\_to\_shuffle))
============================================

temp\_shuffled\_pheno = '../Outputs/Fast-Lmm-Inputs/Temporary-Shuffled-Phenotype.txt'
=====================================================================================

for i in range(NUMBER\_OF\_PERMUTATIONS):
=========================================

time\_permut\_0 = time.time()
=============================

shuffle(indices)
================

phenotype\_shuffled = [](#section-1)
====================================

for j in range(len(indices)):
=============================

phenotype\_shuffled.append(phenotype\_to\_shuffle\[2\]\[indices\[j\]\])
=======================================================================

phenotype\_to\_shuffle\[2\] = phenotype\_shuffled
=================================================

phenotype\_to\_shuffle.to\_csv(temp\_shuffled\_pheno, header=False, index=False, sep=' ')
=========================================================================================

tmp\_shuffled\_df = single\_snp(variants\_to\_test, temp\_shuffled\_pheno,
==========================================================================

cache\_file='../Outputs/Fast-Lmm-Cache/Gwas-Permutations-Cache'+str(i)
----------------------------------------------------------------------

cache\_file='../Outputs/Fast-Lmm-Cache/Gwas-Permutations-Cache.npz',
====================================================================

leave\_out\_one\_chrom=False,
=============================

)
=

tmp\_shuffled\_df\['Full ID'\] = tmp\_shuffled\_df\['Chr'\].astype('str') + '\_' + tmp\_shuffled\_df\['ChrPos'\].astype('str') \# \# \# sorting the new df to match the original \# tmp\_shuffled\_df = tmp\_shuffled\_df\[\['Full ID', 'PValue'\]\] \# tmp\_shuffled\_df\['PValue'\].rename('PValueShuffled'+str(i+1)) \# \# \# results\_df = pd.merge(results\_df, tmp\_shuffled\_df, on='Full ID') \# print('Time for permutation GWAS:' + str(time.time() - time\_permut\_0) + 's')
=======================================================================================================================================================================================================================================================================================================================================================================================================================================================================================

<!-- print 'Time: ' + str(time.clock()-start) -->
<!-- # In[ ]: -->
<!-- # Shuffling ALLELES by VARIANT -->
<!-- for i in range(NUMBER_OF_PERMUTATIONS): -->
<!--     time_permut_0 = time.time() -->
<!--     # Python works a little different than R: Shuffle directly modifies the input data frame! -->
<!--     np.random.shuffle(mysnpdata.val) -->
<!--     Bed.write('VariantsPermuted', mysnpdata) -->
<!--     copyfile(variants_to_test + '.bim', 'VariantsPermuted.bim') -->
<!--     tmp_shuffled_df = single_snp('VariantsPermuted',  phenotype_data, -->
<!-- #                                cache_file='../Outputs/Fast-Lmm-Cache/Gwas-Permutations-Cache'+str(i) -->
<!--                                  cache_file='../Outputs/Fast-Lmm-Cache/Gwas-Permutations-Cache.npz', -->
<!--                                  leave_out_one_chrom=False, -->
<!--                                  ) -->
<!--     tmp_shuffled_df['Full ID'] = tmp_shuffled_df['Chr'].astype('str') + '_' + tmp_shuffled_df['ChrPos'].astype('str') -->
<!--     # sorting the new df to match the original -->
<!--     tmp_shuffled_df = tmp_shuffled_df[['Full ID', 'SNP', 'PValue']] -->
<!--     tmp_shuffled_df = tmp_shuffled_df.rename(columns={'Full ID':'Full IDShuffled'+str(i+1), -->
<!--                                                       'PValue':'PValueShuffled'+str(i+1)}) -->
<!--     snpdata = mysnpdata.val -->
<!--     diff = range(snpdata.shape[1]) -->
<!--     maf = range(snpdata.shape[1]) -->
<!--     n_alleles = range(snpdata.shape[1]) -->
<!--     mean_major = range(snpdata.shape[1]) -->
<!--     for k in range(snpdata.shape[1]): -->
<!--         ref = [j for j, x in enumerate(snpdata[:,k]) if x == 2] -->
<!--         alt = [j for j, x in enumerate(snpdata[:,k]) if x == 0] -->
<!--         meanref = np.mean(pheno[ref]) -->
<!--         meanalt = np.mean(pheno[alt]) -->
<!--         if len(ref) > len(alt): -->
<!--             diff[k] = meanref - meanalt -->
<!--             maf[k] = float(len(alt)) / (len(ref) + len(alt)) -->
<!--             n_alleles[k] = len(ref) + len(alt) -->
<!--             mean_major[k] = meanref -->
<!--         elif len(ref) + len(alt) == 0: -->
<!--             diff[k] = float('NaN') -->
<!--             maf[k] = float('NaN') -->
<!--             n_alleles[k] = len(ref) + len(alt) -->
<!--             mean_major[k] = float('NaN') -->
<!--         else: -->
<!--             diff[k] = meanalt - meanref -->
<!--             maf[k] = float(len(ref)) / (len(ref) + len(alt)) -->
<!--             n_alleles[k] = len(ref) + len(alt) -->
<!--             mean_major[k] = meanalt -->
<!--     diff_df = diff_df = pd.DataFrame(data={'MajMinDiffShuffled'+str(i+1):diff, -->
<!--                                            'MeanMajorShuffled'+str(i+1): mean_major, -->
<!--                                            'NAllelesShuffled'+str(i+1):n_alleles, -->
<!--                                            'MAFShuffled'+str(i+1):maf}) -->
<!--     diff_df['SNP'] = mysnpdata.sid -->
<!--     tmp_shuffled_df = pd.merge(tmp_shuffled_df, diff_df, on='SNP') -->
<!--     tmp_shuffled_df = tmp_shuffled_df.rename(columns={'SNP':'SNPShuffled'+str(i+1)}) -->
<!-- #    results_df = pd.merge(results_df, tmp_shuffled_df, on='Full ID') -->
<!--     results_df = results_df.join(tmp_shuffled_df) -->
<!--     print('Time for permutation GWAS:' + str(time.time() - time_permut_0) + 's') -->
<!-- print 'Time: ' + str(time.clock()-start) -->
<!-- # In[ ]: -->
<!-- results_df.to_csv(OUTPUT_NAME + '-with-Permutations.txt', sep="\t", index=False) -->
<!-- print 'Time: ' + str(time.clock()-start) -->
<!-- # # Manhattan Plot -->
<!-- # In[ ]: -->
<!-- import pylab -->
<!-- import fastlmm.util.util as flutil -->
<!-- flutil.manhattan_plot(results_df.as_matrix(["Chr", "ChrPos", "PValue"]),pvalue_line=1e-5,xaxis_unit_bp=False, plot_threshold=1) -->
<!-- pylab.savefig('../Figures/GWAS_Manhattan_Plot.png') -->
<!-- pylab.close() -->
<!-- print 'Time: ' + str(time.clock()-start) -->
<!-- # # QQ Plot -->
<!-- # In[ ]: -->
<!-- from fastlmm.util.stats import plotp -->
<!-- plotp.qqplot(results_df["PValue"].values, xlim=[0,5], ylim=[0,5]) -->
<!-- pylab.savefig('../Figures/GWAS_QQPlot.png') -->
<!-- pylab.close() -->
<!-- print 'Time: ' + str(time.clock()-start) -->
<!-- # In[ ]: -->
<!-- get_ipython().magic(u'ls') -->
<!-- print 'Time: ' + str(time.clock()-start) -->
<!-- # # The end \o/ -->
<!-- ``` -->
Output image: ![output](Figures/GWAS_Manhattan_Plot.png)

![output](Figures/GWAS_QQPlot.png)

``` r
toc()
```
