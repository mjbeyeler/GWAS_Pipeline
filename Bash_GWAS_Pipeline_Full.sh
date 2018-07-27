#!/bin/bash

# Find out if we're on WSL or not
set -e
if grep -qE "(Microsoft|WSL)" /proc/version &> /dev/null
then
	echo "Windows 10 Bash"
else
	echo "Anything else"
fi



# CONSTANTS

if [ -z "$pheno" ]
then
	pheno="Mass"
fi

if [ -z "$sex" ]
then
        sex="Female"
fi

if [ -z "$maf" ]
then
        maf=0.05
fi

if [ -z "$perm" ]
then
        perm=0
fi

if [ -z "$reproducibility" ]
then
        reproducibility=TRUE
fi

if [ -z "$use_official_gsm" ]
then
        use_official_gsm=FALSE
fi


# Actual Pipeline

Rscript Helper_Scripts/rmd_to_r.R PipelinePart1_AdjustingPhenotypes_BuildingAlleleFilteringScript.Rmd

Rscript PipelinePart1_AdjustingPhenotypes_BuildingAlleleFilteringScript.R $pheno $sex $maf $perm $reproducibility $use_official_gsm

./PipelinePart2_Plink2FilteringAlleles.sh

jupyter nbconvert --to script PipelinePart3_GWASWithPermutations.ipynb
python PipelinePart3_GWASWithPermutations.py 
