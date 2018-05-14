import os
import csv
os.chdir('C:/Users/micha/Documents/Masters-Project-Local/Fast-Lmm-Master')

from fastlmm.association import single_snp

test_snps = "../Data-Import/dgrp2"
test_snps_published = '../Outputs/Plinkfiles/Dgrp2-Tergite-Online'
test_snps_maf5 = '../Outputs/Plinkfiles/dgrp2_maf5'
test_snps_maf5_head = '../Outputs/Plinkfiles/dgrp2_maf5_head'
phenotypes = "../Outputs/Fast-Lmm-Inputs/T5-Pheno-for-Fast-Lmm.txt"
phenotypes_online = '../Outputs/Fast-Lmm-Inputs/T5-Pheno-Online-for-Fast-Lmm.txt'
#cov_fn = "../../tests/datasets/synth/cov.txt"


# Running GWAS
###############################################################################
results_df = single_snp(test_snps_published,  phenotypes_online)

import pylab
import fastlmm.util.util as flutil
flutil.manhattan_plot(results_df.as_matrix(["Chr", "ChrPos", "PValue"]),pvalue_line=1e-5,xaxis_unit_bp=False)
pylab.show()

from fastlmm.util.stats import plotp
plotp.qqplot(results_df["PValue"].values, xlim=[0,5], ylim=[0,5])

# print head of results data frame
import pandas as pd
pd.set_option('display.width', 1000)
results_df.head(n=10)