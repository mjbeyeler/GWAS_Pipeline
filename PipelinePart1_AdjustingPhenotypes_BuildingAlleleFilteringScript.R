## ----loading packages and reproducibility, echo=FALSE, message=FALSE-----

.LIST.OF.PACKAGES <- c(
  'data.table',           #
  'icesTAF',              # dos2unix function
#  'lintr',                # good debugging tool
  'lme4',                 # 
  'reticulate',           # required to switch between R and python chunks, apparently
  'tictoc'                # 
)


source('Helper_Scripts/Environment_Manipulation_and_Reproducibility.R')

# In order to make the script 100% reproducible, keep the next line active:
.BuildReproducibleEnvironment(PROJECT.SNAPSHOT.DATE = "2018-05-16",
                              PROJECT.VERSION       = "3.4.3",
                              SCAN.FOR.PACKAGES     = TRUE)

.LoadPackages(.LIST.OF.PACKAGES)


## ----session info--------------------------------------------------------
cat("The current checkpoint is:\n")
checkpoint::setSnapshot()
cat("\n\nSession Info:\n\n")
sessionInfo()

## ------------------------------------------------------------------------

SEXUAL.DIMORPHISM            <- TRUE
PHENOTYPE.NAME               <- 'Food_Intake'
SEX                          <- 'Female'
NORMALITY.SIGNIFICANCE.LEVEL <- 0.05
INVERSIONS.CONSIDERED        <- c('In.2L.t', 'In.2R.NS', 'In.3R.P', 'In.3R.K', 'In.3R.Mo')

MAF.THRESHOLD=0.05


# GWAS Constants

NUMBER_OF_PERMUTATIONS       <-  0
OUTPUT_NAME                  <-  paste('../Outputs/GWAS', PHENOTYPE.NAME, sep='_')
PHENOTYPE_DATA               <-  paste('../Outputs/Fast-Lmm-Input-', PHENOTYPE.NAME, '-', SEX, '.txt', sep='')
VARIANTS_TO_TEST             <-  '../Outputs/Current_Pipeline_Variants'

fwrite(list(NUMBER_OF_PERMUTATIONS, OUTPUT_NAME, PHENOTYPE_DATA, VARIANTS_TO_TEST),
       file='Outputs/GWAS_Constants', sep="\n")

## ------------------------------------------------------------------------
# Main data

Phenotype_Raw <- read.delim(paste('Inputs/',PHENOTYPE.NAME,'_Phenotype_Full.txt', sep=""), header=T)

# Supporting data

Dgrp2_Inversions <- read.csv('Data/inversion.csv', header=T)
Dgrp2_Infection <- read.csv('Data/wolbachia.csv')

# Some phenotypes actually use the flystock ID instead of the DGRP ID.
# For these cases, this data frame will come in handy.
# Dgrp_Flystock_Ids <- read.delim('Data/Dgrp-Flystocks-Ids.txt',
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
./plink2 --bfile ../Data/dgrp2 --keep ../Outputs/Plink-Lines-$PHENOTYPE_NAME.txt --maf $MAF --make-bed --out ../Outputs/Current_Pipeline_Variants
", sep=''),
file='PipelinePart2_Plink2FilteringAlleles.sh')

# The following command is necessary to 
dos2unix('PipelinePart2_Plink2FilteringAlleles.sh')

