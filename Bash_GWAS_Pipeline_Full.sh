#!/bin/bash

# Find out if we're on WSL or not
set -e
if grep -qE "(Microsoft|WSL)" /proc/version &> /dev/null ; then
    echo "Windows 10 Bash"
else
    echo "Anything else"
fi


Rscript Pipeline_R_Part.R

./Scripts/Plink2_Filtering_Alleles.sh

jupyter nbconvert --to script GWAS_with_Permutations.ipynb
python GWAS_with_Permutations.py 
