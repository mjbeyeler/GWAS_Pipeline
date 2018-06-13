## ------------------------------------------------------------------------
.LIST.OF.PACKAGES <- c(
  'data.table',           #
  'icesTAF',              # dos2unix function
#  'lintr',                # good debugging tool
  'lme4',                 # 
  'reticulate',           # required to switch between R and python chunks, apparently
  'tictoc'                # 
)



INVERSIONS.CONSIDERED        <- c('In.2L.t', 'In.2R.NS', 'In.3R.P', 'In.3R.K', 'In.3R.Mo')
NORMALITY.SIGNIFICANCE.LEVEL <- 0.05
SEXUAL.DIMORPHISM            <- TRUE


if (exists(commandArgs(trailingOnly=T)[1])) {
  .RUN_REPRODUCIBLE            <- commandArgs(trailingOnly=T)[5]
  
  PHENOTYPE.NAME               <- commandArgs(trailingOnly=T)[1]
  SEX                          <- commandArgs(trailingOnly=T)[2]
  
  MAF.THRESHOLD                <- commandArgs(trailingOnly=T)[3]
  
  
  # GWAS Constants
  
  NUMBER_OF_PERMUTATIONS       <- commandArgs(trailingOnly=T)[4]
  USE_OFFICIAL_GSM             <- commandArgs(trailingOnly=T)[6]
} else {
  .RUN_REPRODUCIBLE            <- FALSE
  
  PHENOTYPE.NAME               <- "Mass"
  SEX                          <- "Male"
  
  MAF.THRESHOLD                <- 0.05
  
  
  # GWAS Constants
  
  NUMBER_OF_PERMUTATIONS       <- 0
  USE_OFFICIAL_GSM             <- TRUE
}

OUTPUT_NAME                  <- paste('GWAS', PHENOTYPE.NAME, SEX, sep='_')
PHENOTYPE_DATA               <- paste('Inputs/Fast-Lmm-Input-', PHENOTYPE.NAME, '-', SEX, '.txt', sep='')
VARIANTS_TO_TEST             <- 'Inputs/Current_Pipeline_Variants'

## ----loading packages and reproducibility, echo=FALSE, message=FALSE-----


source('Helper_Scripts/Environment_Manipulation_and_Reproducibility.R')


# In order to make the script 100% reproducible, run_reproducible has to be set to TRUE
run_reproducible <- as.logical(.RUN_REPRODUCIBLE)
if (run_reproducible == TRUE)
  .build_reproducible_environment(PROJECT.SNAPSHOT.DATE = "2018-05-16",
                                  PROJECT.VERSION       = "3.4.3",
                                  SCAN.FOR.PACKAGES     = TRUE)

.load_packages(.LIST.OF.PACKAGES)


## ----session info--------------------------------------------------------
cat("The current checkpoint is:\n")
checkpoint::setSnapshot()
cat("\n\nSession Info:\n\n")
sessionInfo()

## ------------------------------------------------------------------------
# Main data

Phenotype_Raw <- read.delim(paste('Inputs/',PHENOTYPE.NAME,'_Phenotype_Full.txt', sep=""), header=T)

# Supporting data

Dgrp2_Inversions <- read.csv('Raw_Data/inversion.csv', header=T)
Dgrp2_Infection <- read.csv('Raw_Data/wolbachia.csv')

# Some phenotypes actually use the flystock ID instead of the DGRP ID.
# For these cases, this data frame will come in handy.
# Dgrp_Flystock_Ids <- read.delim('Raw_Data/Dgrp-Flystocks-Ids.txt',
#                                 comment.char='#')

## ----functions-----------------------------------------------------------

source('Functions/NormalityHistogram.R')
source('Functions/ChisqForNormality.R')
source('Functions/WriteBare.R')

## ---- echo=FALSE---------------------------------------------------------
source("Scripts/Phenotype_Adjustment.R")

## ------------------------------------------------------------------------
cat(paste("#!/bin/bash

PHENOTYPE_NAME=", PHENOTYPE.NAME,"
MAF=", MAF.THRESHOLD,"

cd plink2_linux_x86_64
./plink2 --bfile ../Raw_Data/dgrp2 --keep ../Inputs/Plink-Lines-$PHENOTYPE_NAME.txt --maf $MAF --make-bed --out ../Inputs/Current_Pipeline_Variants
", sep=''),
file='PipelinePart2_Plink2FilteringAlleles.sh')

# The following command is necessary to 
dos2unix('PipelinePart2_Plink2FilteringAlleles.sh')

## ------------------------------------------------------------------------
fwrite(list(NUMBER_OF_PERMUTATIONS, OUTPUT_NAME, PHENOTYPE_DATA, VARIANTS_TO_TEST, USE_OFFICIAL_GSM),
       file='Inputs/GWAS_Constants.txt', sep="\n")

