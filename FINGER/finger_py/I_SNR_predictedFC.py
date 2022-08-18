# ## All original code is copyright of Graham King, Trinity College Institute of Neuroscience. 
# # Date: 05/12/2021

# Code using predicted SNR-HeadPosition values for two head positions per session (vol 576, vol 1726), 
# in combination with parameters saved from StatsModels Ordinary Least Square linear regressor equations,
# to calculate predicted functional connectivity values.

# Code to calculate the Spearman second order correlation 
# Code to calculate fingerprint identification success

### See the TWO SWITCHES below!!!!!

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



def loadarray():
    snr_all = np.load('/dhcp/fmri_anna_graham/GKgit/snr_npy/snr_96_indivspace.npy')
    mz=snr_all.shape
    print('starting snr_all shape is:')
    print(mz) #snr_all has shape [48, 2, 387]
    print()
    nsubj=mz[0]
    nsess=mz[1]
    nroi=mz[2]
    print ('The dimensions of snr_all is:')
    print(snr_all.ndim)
    print('The size of snr_all is:')
    print(snr_all.size)
    print()

    snr_all_reshaped=np.reshape(snr_all, (nsubj*nsess, nroi))
    print('The shape of snr_all_reshaped is:')
    print(snr_all_reshaped.shape) # snr_all_reshaped has shape [96, 387]
    print()

    # ## Products/edges of SNR values (to match the format of the original snr products used to deduce the OLS parameters c and m):
    # snrbysnr=np.zeros((nsubj*nsess*nhp, nroi, nroi))
    # for ind in range(nsubj*nsess*nhp):
    #     snr=snr_all_reshaped[ind,:]
    #     snrbysnr[ind,:,:]=np.outer(snr, snr) #instead of using snr*snr.T
    # print('The shape of snrbysnr is:')
    # print(snrbysnr.shape) # This is [192,387,387]
    # print()
    # iu1 = np.triu_indices(nroi, k=1)
    # allsnr_iu1=np.array([x[iu1]for x in snrbysnr]) #for each session, selecting the upper triangle of the 387x387 matrix
    # print('The shape of allsnr_iu1 is:')
    # print(allsnr_iu1.shape)
    # print()

    # nedges=allsnr_iu1.shape[1]
    # print('The number of edges is:')
    # print(nedges)
    # print()


    allfc = np.load('/dhcp/fmri_anna_graham/GKgit/finger_npy/fc_96.npy')
    sz=allfc.shape
    print('starting allfc shape is:')
    print(sz)
    print()
    
    allfcmean = allfc.mean(axis=3)
    print('The shape of allfcmean is:')
    print(allfcmean.shape)
    print()
    print('Here are some values of allfcmean:')
    print(allfcmean[:2,:,:5])
    print()
    
    allfcmean_reshaped=np.reshape(allfcmean, (nsubj*nsess, nroi))
    print('The shape of allfcmean_reshaped is:')  # This shape is (96, 387)
    print(allfcmean_reshaped.shape)
    print()

    return snr_all_reshaped, allfcmean_reshaped, nroi


def predicted_fc(snr_all_reshaped, nroi):

    parameters=np.load('/dhcp/fmri_anna_graham/GKgit/snr_npy/params_416_indivspace.npy') #Edit name here
    print('The shape of parameters array is:')
    print(parameters.shape) #This has shape (387, 2)
    print()

    predicted=np.zeros((96, nroi))
    for parind, (c, m) in enumerate(parameters):
        ## y = c + mx
        fit = c + m*snr_all_reshaped[:, parind]
        predicted[:, parind]=fit
    print('The shape of predicted array is:')
    print(predicted.shape) # should be (96, 387)
    print()
    print('Some values of the predicted array are:')
    print(predicted[0:2,0:10:1])
    print()

    return predicted


def compare(predicted, allfcmean_reshaped):
    
    fc_corr=np.corrcoef(allfcmean_reshaped, predicted, rowvar = True)
    print('the shape of the fc_corr matrix is:')
    print(fc_corr.shape)
    print()

    plt.figure(1)
    plt.imshow(fc_corr)
    plt.colorbar()
    plt.savefig('/dhcp/fmri_anna_graham/GKgit/head_position/FINGER/finger_figures/fingerfigures_misc/fc_corr.jpg')

    compareFC=fc_corr[0:96,96:192]
    print('The shape of compareFC is:')
    print(compareFC.shape)
    print()

    plt.figure(2)
    plt.imshow(compareFC)
    plt.colorbar()
    plt.savefig('/dhcp/fmri_anna_graham/GKgit/head_position/FINGER/finger_figures/fingerfigures_misc/FCreal_FCpredicted.jpg')


if __name__ == '__main__':

    snr_all_reshaped, allfcmean_reshaped, nroi = loadarray()
    predicted = predicted_fc(snr_all_reshaped, nroi)
    compare(predicted, allfcmean_reshaped)

