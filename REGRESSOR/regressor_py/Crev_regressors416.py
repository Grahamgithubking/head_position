#### All original code is copyright of Graham King, Trinity College Institute of Neuroscience.
### Using predicted SNR-HeadPosition and corresponding Fuctional Connectivity data to calculate the parameters of 74,691 linear regressors across 416 participants with one session.
### Note: The raw SNR data from these 416 participants was used to create the SNR-Coil (SNR of the Coil)

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


def Loading_reshaping(snrcoil, snrtrue):
    #Loading and Reshaping both predicted SNR-HeadPosition data and FC data
    
    allfc = np.load('/dhcp/fmri_anna_graham/GKgit/finger_npy/fc_416.npy')
    sz=allfc.shape
    print('starting allfc shape is:')
    print(sz)
    nsubj=sz[0]
    nsess=sz[1]
    nroi=sz[2]
    
    iu1=np.triu_indices(nroi, k=1) #selecting the upper triangle of a matrix

    allfc_reshaped=np.reshape(allfc, (nsubj*nsess, nroi, nroi))
    print('The shape of allfc_reshaped is:')
    print(allfc_reshaped.shape)
    allfc_iu1=np.array([x[iu1]for x in allfc_reshaped]) #for each session, selecting the upper triangle of the 387x387 matrix
    print('The shape of allfc_iu1 is:')
    print(allfc_iu1.shape)
    print('Some values of allfc_iu1 are:')
    print(allfc_iu1[:2,:5])
    print()

    # Calculating the number of edges:
    nedges=allfc_iu1.shape[1]
    print('The number of edges is:')
    print(nedges)
    print()

    allfc = allfc_iu1

    # # Getting the Fishers r2z standardization:
    # allfc_fz=np.arctanh(allfc_iu1)
    # print('allfc_fz has shape:')
    # print(allfc_fz.shape)
    # print('Some values of allfc_fz are:')
    # print(allfc_fz[:2,:5])
    # print()

    ######################################

    if snrcoil:
        snr_all = np.load('/dhcp/fmri_anna_graham/GKgit/snr_npy/416_snr_416_erode1.npy')
    if snrtrue:
        snr_all = np.load('/dhcp/fmri_anna_graham/GKgit/snr_npy/416_snr_true.npy')

    mz=snr_all.shape
    nsnr=mz[2]

    snr_all_reshaped=np.reshape(snr_all, (nsubj*nsess, nsnr))
    print('The shape of snr_all_reshaped is:')
    print(snr_all_reshaped.shape)
    print()

    # Multiplying the SNR across all edges for each ROI:
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
    print('Some values of allsnr_iu1 are:')
    print(allsnr_iu1[:2,:5])
    print()

    allsnr = allsnr_iu1

    # # Standardizing the SNR edge values:
    # smean=np.mean(allsnr_iu1, axis=1)
    # sstd=np.std(allsnr_iu1, axis=1)
    # allsnr_sz=np.zeros((nsubj*nsess,nedges))
    # for ind in range(nsubj*nsess):
    #     allsnr_sz[ind,:]=(allsnr_iu1[ind,:]-smean[ind,])/sstd[ind,]
    # print('Some values of allsnr_sz are:')
    # print(allsnr_sz[:2,:5])
    # print()


    return allfc, allsnr, nedges

def Edge_regressors(allfc, allsnr, nedges):
    ## Regression models for EACH edge (across 416 participants with < 1 > session)
     
    parameters=np.zeros((nedges,2)) #nedges OR specify number nedge of 74,691 here!
    p=np.zeros(nedges)
    r2=np.zeros(nedges)

    for r in range(nedges):
        x=allsnr[:,r]
        y=allfc[:,r]

        X2 = sm.add_constant(x)
        results = sm.OLS(y, X2).fit()
        ## y = c + mx where (c, m) are parameters
        ## NOTE: param1 is the constant_value, param2 is the coefficient_x1
        parameters[r] = results.params
        p[r] = results.pvalues[1]
        r2[r] = results.rsquared

    return parameters, p, r2

    
    ## NOTE: param1 is the constant_value, param2 is the coefficient_x1
    
    
if __name__ == '__main__':

    snrcoil = True
    snrtrue = False

    allfc, allsnr, nedges = Loading_reshaping(snrcoil, snrtrue)

    parameters, p, r2 = Edge_regressors(allfc, allsnr, nedges)

    if snrcoil:
        np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/params_416_erode1.npy', parameters)
        np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/p_416_erode1.npy', p)
        np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/r2_416_erode1.npy', r2)
    if snrtrue:
        np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/params_416_snrtrue.npy', parameters)

    
    
