import os
#import csv
os.chdir('C:/Users/micha/OneDrive/Desktop/Master\'s-Project/FaST-LMM-master')

from fastlmm.association import single_snp

test_snps = "../Data/dgrp2"
test_snps_maf5 = '../Outputs/Dgrp2-Plinkfiles/dgrp2_maf5'
test_snps_maf5_head = '../Outputs/Dgrp2-Plinkfiles/dgrp2_maf5_first1000'
phenotype_t5 = '../Outputs/Fast-Lmm-Inputs//Fast-Lmm-Input-Tergites-Female.txt'
phenotype_female = "../Outputs/Fast-Lmm-Inputs/Fast-Lmm-Input-Mass-Female.txt"
phenotype_male = "../Outputs/Fast-Lmm-Inputs/Fast-Lmm-Input-Mass-Male.txt"
phenotype_dimorphism = "../Outputs/Fast-Lmm-Inputs/Fast-Lmm-Input-Mass-Dimorphism.txt"
phenotypes_online = '../Outputs/Fast-Lmm-Inputs/T5-Pheno-Online-for-Fast-Lmm.txt'
#cov_fn = "../../tests/datasets/synth/cov.txt"


# Running GWAS
###############################################################################
results_df_female = single_snp(test_snps_maf5,  phenotype_dimorphism,
#                               leave_out_one_chrom=False,
                               output_file_name='../Outputs/Fast-Lmm-Outputs/Gwas-Maf5-Mass-Dimorphism.txt')
#results_df_male = single_snp(test_snps_maf5,  phenotype_male, leave_out_one_chrom=False)

import pylab
import fastlmm.util.util as flutil
flutil.manhattan_plot(results_df_female.as_matrix(["Chr", "ChrPos", "PValue"]),pvalue_line=1e-5,xaxis_unit_bp=False)
pylab.show()

from fastlmm.util.stats import plotp
plotp.qqplot(results_df_female["PValue"].values, xlim=[0,5], ylim=[0,5])

# print head of results data frame
import pandas as pd
pd.set_option('display.width', 1000)
results_df_female.head(n=10)