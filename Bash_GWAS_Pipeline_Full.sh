#!/bin/bash

# Find out if we're on WSL or not
set -e
if grep -qE "(Microsoft|WSL)" /proc/version &> /dev/null ; then
    echo "Windows 10 Bash"
else
    echo "Anything else"
fi


Rscript Helper_Scripts/Rmd_to_R.R PipelinePart1_AdjustingPhenotypes_BuildingAlleleFilteringScript.R

Rscript PipelinePart1_AdjustingPhenotypes_BuildingAlleleFilteringScript.R

./PipelinePart2_Plink2FilteringAlleles.sh

jupyter nbconvert --to script PipelinePart3_GWASWithPermutations.ipynb
python PipelinePart3_GWASWithPermutations.py 
