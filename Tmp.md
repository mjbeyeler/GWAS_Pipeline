GWAS Pipeline Full
================
Michael Beyeler
2018-05-13

Setup
=====

``` r
source('Helper_Scripts/Reprobucibility_Functions.R')
.BuildReproducibleEnvironment('2018-01-01')
```

    ## Scanning for packages used in this project

    ## No file at path 'C:\Users\micha\AppData\Local\Temp\RtmpC0ufeP\file363c7bd92861.Rmd'.

    ## No file at path 'C:\Users\micha\AppData\Local\Temp\RtmpC0ufeP\file363c37443c1b.Rmd'.

    ## - Discovered 5 packages

    ## Unable to parse 2 files:

    ## - GWAS_Pipeline_Full.Rmd

    ## - Tmp.Rmd

    ## Packages not available in repository and won't be installed:

    ##  - not.a.package

    ##  - RevoUtils

    ##  - RevoUtilsMath

    ## All detected packages already installed

    ## checkpoint process complete

    ## ---

``` r
.LIST.OF.PACKAGES <- c(
  'data.table',           #
  'tictoc',               # 
  'lme4',                 # 
  'reticulate',           # used to properly being able to switch between python and r in R Markdown
  'icesTAF'               # dos2unix function
)
.LoadPackages(.LIST.OF.PACKAGES)
```

    ## $data.table
    ## NULL
    ## 
    ## $tictoc
    ## NULL
    ## 
    ## $lme4
    ## NULL
    ## 
    ## $reticulate
    ## NULL
    ## 
    ## $icesTAF
    ## NULL

``` r
# Functions

source('Functions/NormalityHistogram.R')
source('Functions/ChisqForNormality.R')
source('Functions/WriteBare.R')


# CUSTOM SESSION PYTHON PATH, MUST BE ANACONDA2 PYTHON
knitr::opts_chunk$set(engine.path = list(python = '/anaconda/bin/python'))
knitr::knit_engines$set(python = '/anaconda/bin/python')

sessionInfo()
```

    ## R version 3.4.3 (2017-11-30)
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
    ## [1] icesTAF_1.4-1       reticulate_1.3.1    lme4_1.1-15        
    ## [4] Matrix_1.2-12       tictoc_1.0          data.table_1.10.4-3
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_0.12.14     knitr_1.18       magrittr_1.5     splines_3.4.3   
    ##  [5] MASS_7.3-47      lattice_0.20-35  R6_2.2.2         minqa_1.2.4     
    ##  [9] httr_1.3.1       stringr_1.2.0    tools_3.4.3      grid_3.4.3      
    ## [13] nlme_3.1-131     htmltools_0.3.6  yaml_2.1.16      checkpoint_0.4.3
    ## [17] rprojroot_1.3-1  digest_0.6.13    nloptr_1.0.4     evaluate_0.10.1 
    ## [21] rmarkdown_1.6    stringi_1.1.6    compiler_3.4.3   backports_1.1.2 
    ## [25] jsonlite_1.5

FaST-LMM GWAS
=============

``` python
import time
start = time.clock()

import os
os.chdir('FaST_LMM')

from fastlmm.association import single_snp
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

print str(time.clock-start)
```
