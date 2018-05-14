#!/bin/bash

PHENOTYPE_NAME='Mass'
MAF=0.05

cd plink2_linux_x86_64
./plink2 --bfile ../Data/dgrp2 --keep ../Outputs/Plink-Lines-$PHENOTYPE_NAME.txt --maf $MAF --make-bed --out ../Outputs/MassVariants_MAF5
