#### All original code is copyright of Graham King, Trinity College Institute of Neuroscience.

### Using predicted SNR-HeadPosition and corresponding Fuctional Connectivity data to calculate the parameters of 74,691 linear regressors across 48 participants with two sessions
### Note: The raw SNR data from these 48 participants was NOT used to create the SNR-Coil (SNR of the Coil). This is an independant group.

import pandas as pd
import numpy as np
import glob
import os
import nilearn.signal
import matplotlib.pyplot as plt
from scipy import stats
import pingouin as pg
import statsmodels.api as sm
import seaborn as sns


def Loading_reshaping():
    
    allfc = np.load('/dhcp/fmri_anna_graham/GKgit/finger_npy/fc_96.npy')
    sz=allfc.shape
    print('starting allfc shape is:')
    print(sz)
    print()
    nsubj=sz[0]
    nsess=sz[1]
    nroi=sz[2]
    
    iu1=np.triu_indices(nroi, k=1) #for each session, selecting the upper triangle of the 387x387 matrix

    allfc_reshaped=np.reshape(allfc, (nsubj*nsess, nroi, nroi))
    print('The shape of allfc_reshaped is:')
    print(allfc_reshaped.shape)
    allfc_iu1=np.array([x[iu1]for x in allfc_reshaped]) #for each session, selecting the upper triangle of the 387x387 matrix
    print('The shape of allfc_iu1 is:')
    print(allfc_iu1.shape)
    print()

    snr_all = np.load('/dhcp/fmri_anna_graham/GKgit/snr_npy/snr_96_erode.npy')
    mz=snr_all.shape
    nsnr=mz[2]
    snr_all_reshaped=np.reshape(snr_all, (nsubj*nsess, nsnr))
    print('The shape of snr_all_reshaped is:')
    print(snr_all_reshaped.shape)
    print()
    snrbysnr=np.zeros((nsubj*nsess,nsnr,nsnr))
    for ind in range(nsubj*nsess):
        snr=snr_all_reshaped[ind,:]
        snrbysnr[ind,:,:]=np.outer(snr, snr) #instead of using snr*snr.T
    print('The shape of snrbysnr is:')
    print(snrbysnr.shape)
    print()
    allsnr_iu1=np.array([x[iu1]for x in snrbysnr]) #for each session, selecting the upper triangle of the 387x387 matrix
    print('The shape of allsnr_iu1 is:')
    print(allsnr_iu1.shape)

    nedges=allfc_iu1.shape[1]
    print('The number of edges is:')
    print(nedges)

    return allfc_iu1, allsnr_iu1, nedges

def Edge_regressors(allfc_iu1, allsnr_iu1, nedges):
    ## Regression models for EACH edge (across 48 participants with < 2 > sessions):

    parameters=np.zeros((nedges,2)) #nedges
    p=np.zeros(nedges)
    r2=np.zeros(nedges)

    for r in range(nedges): #nedges OR specify number nedge of 74,691 here!
        x=allsnr_iu1[:,r]
        y=allfc_iu1[:,r]

        X2 = sm.add_constant(x)
        results = sm.OLS(y, X2).fit()
        ## y = c + mx where (c, m) are parameters
        ## NOTE: param1 is the constant_value, param2 is the coefficient_x1
        parameters[r] = results.params
        p[r] = results.pvalues[1]
        r2[r] = results.rsquared


    print('The parameters array size is:')
    print(parameters.shape)
    print()
    ## NOTE: param1 is the constant_value, param2 is the coefficient_x1
    np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/params_96_erode.npy', parameters)


    
if __name__ == '__main__':

    allfc_iu1, allsnr_iu1, nedges = Loading_reshaping()

    Edge_regressors(allfc_iu1, allsnr_iu1, nedges)

    



