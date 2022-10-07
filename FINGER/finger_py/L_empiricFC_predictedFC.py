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

    snr_all = np.load('/dhcp/fmri_anna_graham/GKgit/snr_npy/96_snr_416_erode1.npy')
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

def predicted_fc(allsnr_iu1, nedges):
    ### Calculating the predicted fc using OLS parameters and edge-level products of predicted SNR-group values

    parameters=np.load('/dhcp/fmri_anna_graham/GKgit/snr_npy/74691_params_snr_416_erode1.npy') # parameters from 416 sessions
    print('The shape of parameters array is:')
    print(parameters.shape)
    print()

    predicted=np.zeros((96, nedges))
    for parind, (c, m) in enumerate(parameters):
        ## y = c + mx
        fit = c + m*allsnr_iu1[:, parind]
        predicted[:, parind]=fit
    ## predicted array size [96, 74691]
    print('Some values of predicted array are:')
    print(predicted[0:2,0:10])
    print()

    return predicted

def correlation(predicted, allfc_iu1):

    ## Second order Spearman correlation:
    correl = np.corrcoef(predicted, allfc_iu1)
    print('The shape of the correl matrix is:')
    print(correl.shape)
    print()

    compare=correl[0:96:1,96::1] # The upper right quadrant of the 192x192 matrix
    diagonal=np.diag(compare) # The diagonal = This is the correlation values between empirical vs predicted FC for each of 96 sessions.
    print(f"The shape of diagonal is: {diagonal.shape}")
    print(f"Some values of diagonal are:")
    print(diagonal[:10])
    print()

    print('The mean of the diagonal values is:')
    print(np.mean(diagonal))
    print()

    plt.figure(1)
    sns.distplot(diagonal, kde=True) 
    plt.xlabel('Correlation (Pearson)')
    plt.savefig('/dhcp/fmri_anna_graham/GKgit/head_position/FINGER/finger_figures/fingerfigures_misc/emp_v_pred.jpg') # Edit name here!



if __name__ == '__main__':

    allfc_iu1, allsnr_iu1, nedges = Loading_reshaping()

    predicted = predicted_fc(allsnr_iu1, nedges)

    correlation(predicted, allfc_iu1)
