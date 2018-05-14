from fastlmm.inference.lmm_cov import LMM
import numpy as np

Gsm = np.genfromtxt('../Data/freeze2.common.rel.mat', skip_header=True)
Gsm = Gsm[:, 2:207]
Gsm = LMM(K=Gsm)
Gsm.setSU_fromK()

np.savetxt('../Outputs/Fast-Lmm-Inputs/Ncsu-Gsm-U.txt', Gsm.U)
np.savetxt('../Outputs/Fast-Lmm-Inputs/Ncsu-Gsm-S.txt', Gsm.S)