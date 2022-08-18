### to be continued....

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


def loadsnr(imageroot, one_sessid, two_sessid, snrcoil, snrtrue, split):

    if one_sessid:
        x='snrcoil'
        y='416'
        snr_all=np.load('/dhcp/fmri_anna_graham/GKgit/snr_npy/416_snr_416_erode1.npy') #This has shape [416,1,387]

        sz=snr_all.shape
        print('The shape of snr_all is:')
        print(sz)
        print()
        nsubj=sz[0]
        nsess=sz[1]
        nroi=sz[2]

        snr_all_reshaped=np.reshape(snr_all, (nsubj*nsess, nroi))
        print('The shape of snr_all_reshaped is:')
        print(snr_all_reshaped.shape)
        print()

    if two_sessid:
        if snrcoil:
            x='snrcoil'
            y='96'
            snr_all=np.load('/dhcp/fmri_anna_graham/GKgit/snr_npy/96_snr_416_erode1.npy') #This has shape [48,2,387]
        if snrtrue:
            x='snrtrue'
            y='96'
            snr_all=np.load('/dhcp/fmri_anna_graham/GKgit/snr_npy/96_snr_true.npy') #This has shape [48,2,387]

        sz=snr_all.shape
        print('The shape of snr_all is:')
        print(sz)
        print()
        nsubj=sz[0]
        nsess=sz[1]
        nroi=sz[2]
        
        snr_all_reshaped=np.reshape(snr_all, (nsubj*nsess, nroi))
        print('The shape of snr_all_reshaped is:')
        print(snr_all_reshaped.shape)
        print()

    if split:
        if snrcoil:
            x='snrcoil'
            y='192'
            snr_all = np.load('/dhcp/fmri_anna_graham/GKgit/snr_npy/192hp_snr.npy')
        if snrtrue:
            x='snrtrue'
            y='192'
            snr_all = np.load('/dhcp/fmri_anna_graham/GKgit/snr_npy/192hp_snrtrue.npy')        

        sz=snr_all.shape
        print('The shape of snr_all is:')
        print(sz)
        print()
        nsubj=sz[0]
        nsess=sz[1]
        nhp=sz[2]
        nroi=sz[3]

        snr_all_reshaped=np.reshape(snr_all, (nsubj*nsess*nhp, nroi))
        print('The shape of snr_all_reshaped is:')
        print(snr_all_reshaped.shape)
        print()

    # I <did/did not> DO THIS:
    ##### Multiplying the SNR across all edges for each ROI:
    # snrbysnr=np.zeros((nsubj*nsess,nroi,nroi))
    # for ind in range(nsubj*nsess):
    #     snr=snr_all_reshaped[ind,:]
    #     snrbysnr[ind,:,:]=np.outer(snr, snr) #instead of using snr*snr.T
    # print('The shape of snrbysnr is:')
    # print(snrbysnr.shape)
    # print()
    # iu1=np.triu_indices(nroi, k=1) #selecting the upper triangle of a matrix
    # allsnr_iu1=np.array([x[iu1]for x in snrbysnr]) #for each session, selecting the upper triangle of the 387x387 matrix
    # print('The shape of allsnr_iu1 is:')
    # print(allsnr_iu1.shape)
    # snr_all_reshaped=allsnr_iu1


    snrim, pval=stats.spearmanr(snr_all_reshaped.T)
    print('the shape of the spearman correlation of snr_all is:')
    print(snrim.shape)
    print()           

    plt.figure(1)
    plt.imshow(snrim)
    plt.colorbar()
    plt.title('SNR Correlation_' + str(y) + '_' + str(x))
    plt.savefig(str(imageroot) + str(y) + str(x) + '.jpg')


if __name__ == '__main__':

    imageroot='/dhcp/fmri_anna_graham/GKgit/head_position/FINGER/finger_figures/fingerfigures_misc/'

    ## SWITCHES:
    one_sessid=False
    two_sessid=True

    split=False
    snrcoil=True
    snrtrue=False
  
    loadsnr(imageroot, one_sessid, two_sessid, snrcoil, snrtrue, split)
