import pandas as pd
import numpy as np
import glob
import os
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
from scipy import stats
import seaborn as sns


# Did you use the 7 or the 17 network Schaefer? If the 7 then why not pick these regions (I've chosen region 1 for each, use something else if you can think of a good reason):
# 149: 7Networks_LH_Default_Temp_1,
# 159: 7Networks_LH_Default_Par_1,    
# 166: 7Networks_LH_Default_PFC_1,   
# 190: 7Networks_LH_Default_pCunPCC_1
# and same for RH, then select the FC and SNR and scatter plot like this
# regions = [149, 159, 166, 190, ...]
# fc=fc[regions,:][:,regions] # add subjects : too
# snr=snr[regions,:][:,regions]
# plt.scatter(np.ravel(snr),np.ravel(fc)

# 154	17Networks_LH_DefaultA_pCunPCC_1
# 161	17Networks_LH_DefaultA_PFCm_1
# 167	17Networks_LH_DefaultB_Temp_1
# 175	17Networks_LH_DefaultB_PFCd_1
# 363	17Networks_RH_DefaultA_pCunPCC_1
# 368	17Networks_RH_DefaultA_PFCm_1
# 374	17Networks_RH_DefaultB_Temp_1
# 377	17Networks_RH_DefaultB_PFCd_1

rois = np.array([154,161,167,175,363,368,374,377]) # base 1
rois = rois -1 # base 0
print(f"These are the zero'd regions: {rois}")

### Loading allfc:
fc = np.arange(48)
fc_reshaped=np.reshape(fc, (3,4,4))
print(fc_reshaped)
nsubj = np.arange(3)

### Loading allsnr:
tbc....

def interactions(item_reshaped, regions):
    mat = np.empty([len(nsubj),len(regions),len(regions)])
    for subind, subj in enumerate(nsubj):
        for inda, regiona in enumerate(regions):
            for indb, regionb in enumerate(regions):
                item_temp = item_reshaped[subj,regiona,regionb]
                mat[subind,inda,indb] = item_temp
    print(f"These are some values of mat: {mat[:3,:4,:4]}")
    print()
    iu1=np.triu_indices(len(regions), k=1) #selecting the upper triangle of a matrix (not including main diagonal)
    mat_iu1=np.array([x[iu1]for x in mat])

    return mat_iu1

fc_iu1 = interactions(fc_reshaped, rois)
snr_iu1 = interactions(snr_reshaped, rois)

print(f"The shape of fc_iu1 is: {fc_iu1.shape}") ## This should be shape [416, 8] ???
print()

plt.figure(1)
plt.scatter(fc_iu1, snr_iu1)
plt.savefig('/dhcp/fmri_anna_graham/GKgit/head_position/FINGER/finger_figures/fingerfigures_misc/SNR_covary_FC.jpg')
