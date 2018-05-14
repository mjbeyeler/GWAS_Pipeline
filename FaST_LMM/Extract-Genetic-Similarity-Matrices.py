from fastlmm.association import gsm_extraction_leave_out_one_chrom_false
#from fastlmm.inference.lmm_cov import LMM

# Containing all nuclear variants filtered by Phred quality score 500 and other
# (more details in Huang2014)
variants_nuclear = '../Data/dgrp2'

variants_mito = '../Outputs/Plinkfiles/MitoSeq_AllRuns_dm6_chrM.annot.biallellic_ConvertedReference'

# A randomly generated phenotype, containing all DGRP2 lines.
random_artificial_pheno = '../Outputs/Fast-Lmm-Inputs/Dgrp2-AllLines-RandomPheno.txt'
random_artificial_pheno_mito = '../Outputs/Fast-Lmm-Inputs/Dgrp2-AllLines-RandomPheno-for-Mito.txt'
random_artificial_pheno_mito_with_non_dgrp = '../Outputs/Fast-Lmm-Inputs/Dgrp2-AllLines-RandomPheno-for-Mito-IncludingNonDgrp.txt'
tergite5_pheno = '../Outputs/Fast-Lmm-Inputs/Fast-Lmm-Input-Tergites-Female.txt'



# Running a GWAS with storing cache, in order to output the spectrally
# transformed GSM.
gsm_extraction_leave_out_one_chrom_false(variants_nuclear,
           random_artificial_pheno,
           cache_file = '../Outputs/GSM-NCSUFreeze2Full')
#
#gsm_extraction_leave_out_one_chrom_false(variants_nuclear,
#           tergite5_pheno,
#           cache_file = '../Outputs/Fast-Lmm-Outputs/GsmNuclear-Tergite.txt')


#lmm_nuclear = LMM()
#with np.load('../Outputs/Fast-Lmm-Cache/Cache-Nuclear.npz') as data:
#    nuclear_U = data['arr_0']
#    nuclear_S_vector = data['arr_1']

#nuclear_S = np.diag(nuclear_S_vector)
#
#nuclear_K = np.dot(np.dot(nuclear_U, nuclear_S), nuclear_U)
#            
#nmm_nuclear.get_K


# The following code has to be ran MANUALLY
gsm_extraction_leave_out_one_chrom_false(variants_mito,
           random_artificial_pheno_mito,
           cache_file = '../Outputs/Fast-Lmm-Outputs/GsmMito-171DgrpLines.txt')
gsm_extraction_leave_out_one_chrom_false(variants_mito,
           random_artificial_pheno_mito_with_non_dgrp,
               cache_file = '../Outputs/Fast-Lmm-Outputs/GsmMito-All179Lines.txt')