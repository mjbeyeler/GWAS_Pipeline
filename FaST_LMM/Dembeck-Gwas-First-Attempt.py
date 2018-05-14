import os
import csv
os.chdir('C:/Users/micha/OneDrive/Desktop/Master\'s-Project/Fast-Lmm-Master')

from fastlmm.association import single_snp

test_snps = "../Data-Import/dgrp2"
test_snps_maf5 = '../Data/Dgrp2-Plinkfiles/dgrp2_maf5'
test_snps_maf5_head = '../Data/Dgrp2-Plinkfiles/dgrp2_maf5_head'
phenotypes = "../Data/Fast-Lmm-Inputs/T5-Pheno-for-Fast-Lmm.txt"
phenotypes_online = '../Data/Fast-Lmm-Inputs/T5-Pheno-Online-for-Fast-Lmm.txt'
#cov_fn = "../../tests/datasets/synth/cov.txt"


# Running GWAS with only the first 100 SNPs
###############################################################################
results_first = single_snp(test_snps_maf5_head,  phenotypes)
results_first_online = single_snp(test_snps_maf5_head,  phenotypes_online)

# displaying results
#results_first['PValue']
results_first['PValue']

# writing results to csv
csvfile = "../Data/Fast-Lmm-Outputs/D15-Gwas-First26.csv"
csvfile_online = "../Data/Fast-Lmm-Outputs/D15-Gwas-First26-Online.csv"
# This would work for simpler things
#with open(csvfile, "w") as output:
#    writer = csv.writer(output, lineterminator='\n')
#    [writer.writerow(r) for r in results_first]    
#
##Assuming res is a list of lists
#with open(csvfile, "w") as output:
#    writer = csv.writer(output, lineterminator='\n')
#    writer.writerows(results_first)

# Results are a Pandas dataframe, so I can just use the to_csv feature:
results_first.to_csv(csvfile)
results_first_online.to_csv(csvfile_online)