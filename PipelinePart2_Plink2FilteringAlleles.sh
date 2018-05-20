#!/bin/bash

PHENOTYPE_NAME=Food_Intake
MAF=0.05

cd plink2_linux_x86_64
./plink2 --bfile ../Data/dgrp2 --keep ../Outputs/Plink-Lines-$PHENOTYPE_NAME.txt --maf $MAF --make-bed --out ../Outputs/Current_Pipeline_Variants