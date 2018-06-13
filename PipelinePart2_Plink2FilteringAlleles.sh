#!/bin/bash

PHENOTYPE_NAME=CS
MAF=0.4

cd plink2_linux_x86_64
./plink2 --bfile ../Raw_Data/dgrp2 --allow-extra-chr --keep ../Inputs/Plink-Lines-$PHENOTYPE_NAME.txt --maf $MAF --make-bed --out ../Inputs/Current_Pipeline_Variants
