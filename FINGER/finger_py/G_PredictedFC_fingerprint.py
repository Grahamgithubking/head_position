# ## All original code is copyright of Graham King, Trinity College Institute of Neuroscience. 
# ## Last edited Oct 2022

# Code using predicted SNR-HeadPosition values for two head positions per session (vol 537, vol 1763), 
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
    # Loading predicted SNR-Headposition values for two head positions per session (vol 537, vol 1763),
    # for all sessions of the 44 particpants with both preterm and term sessions
    snr_all = np.load('/dhcp/fmri_anna_graham/GKgit/snr_npy/176hp_snr.npy')
    mz=snr_all.shape
    print('starting snr_all shape is:')
    print(mz) #snr_all has shape [44, 2, 2, 387]
    print()
    nsubj=mz[0]
    nsess=mz[1]
    nhp=mz[2]
    nroi=mz[3]
    print ('The dimensions of snr_all is:')
    print(snr_all.ndim)
    print('The size of snr_all is:')
    print(snr_all.size)
    print()

    snr_all_reshaped=np.reshape(snr_all, (nsubj*nsess*nhp, nroi))
    print('The shape of snr_all_reshaped is:')
    print(snr_all_reshaped.shape) # snr_all_reshaped has shape [176, 387]
    print()

    ## Products/edges of SNR values (to match the format of the original snr products used to deduce the OLS parameters c and m):
    snrbysnr=np.zeros((nsubj*nsess*nhp, nroi, nroi))
    for ind in range(nsubj*nsess*nhp):
        snr=snr_all_reshaped[ind,:]
        snrbysnr[ind,:,:]=np.outer(snr, snr) #instead of using snr*snr.T
    print('The shape of snrbysnr is:')
    print(snrbysnr.shape) # This is [176,387,387]
    print()
    iu1 = np.triu_indices(nroi, k=1)
    allsnr_iu1=np.array([x[iu1]for x in snrbysnr]) #for each session, selecting the upper triangle of the 387x387 matrix
    print('The shape of allsnr_iu1 is:')
    print(allsnr_iu1.shape)
    print()

    nedges=allsnr_iu1.shape[1]
    print('The number of edges is:')
    print(nedges)
    print()

    return allsnr_iu1, nsubj, nsess, nhp, nedges


def predicted_fc(allsnr_iu1, nsubj, nsess, nhp, nedges):
    ### Calculating the predicted fc using OLS parameters and edge-level products of predicted SNR-HeadPosition values for 2 head positions

    parameters=np.load('/dhcp/fmri_anna_graham/GKgit/snr_npy/params_416_erode1.npy') #Edit name here
    print('The shape of parameters array is:')
    print(parameters.shape)
    print()

    predicted=np.zeros((nsubj*nsess*nhp, nedges))
    for parind, (c, m) in enumerate(parameters):
        ## y = c + mx
        fit = c + m*allsnr_iu1[:, parind]
        predicted[:, parind]=fit
    print('The shape of predicted array is:')
    print(predicted.shape) # should be (176, 74691)
    print('Some values are:')
    print(predicted[0,0:10:1])
    print()

    return predicted

def fingerprint(oldersession, youngersession, acrosssession, predicted, nsubj):

    ## Second order Spearman correlation:
    sessim, pval=stats.spearmanr(predicted.T)
    print('the shape of the sessim matrix is:')
    print(sessim.shape)
    print()

    #Find all the points comparing across subjects and sessions for 1 head position, position volume 537 ONLY!:
    if acrosssession:
        comparesess=sessim[2::4,0::4] #this chooses across session values for position 537
        print('the shape of across session comparehp is:')
        print(comparesess.shape)
        print()

        plt.figure(figsize=(12,8))
        c = plt.imshow(comparesess, cmap='Greens') # Edit color here
        plt.colorbar(c)
        x = np.arange(1,45,1) # the grid to which your data corresponds
        nx = x.shape[0]
        no_labels = 11 # how many labels to see on axis x
        step_x = int(nx / (no_labels - 1)) # step between consecutive labels
        x_positions = np.arange(3,nx+1,step_x) # pixel count at label position
        x_labels = x[3::step_x] # labels you want to see
        plt.xticks(x_positions, x_labels)
        plt.yticks(x_positions, x_labels)
        plt.title('RSM of Predicted FC values across sessions \n (using one head position)', fontsize=20, fontweight="bold") # Edit title
        plt.xlabel('Participants No.', fontsize=18)
        plt.ylabel('Participants No.', fontsize=18)
        plt.savefig('/dhcp/fmri_anna_graham/GKgit/head_position/FINGER/finger_figures/predFC_figures/RSM_predFC_acrosssess.png')


    #Find all the points comparing across subjects and across head positions (HP537 vs HP1763) within each session:
    if oldersession:
        comparesess=sessim[2::4,3::4] #this chooses the older session/timepoint
        print('the shape of older comparehp is:')
        print(comparesess.shape)
        print()
    if youngersession:
        comparesess=sessim[0::4,1::4] #this chooses the younger session/timepoint
        print('the shape of younger comparehp is:')
        print(comparesess.shape)
        print()

    if oldersession:
        plt.figure(figsize=(12,8))
        c = plt.imshow(comparesess, cmap='Reds') # Edit color here
        plt.colorbar(c)
        x = np.arange(1,45,1) # the grid to which your data corresponds
        nx = x.shape[0]
        no_labels = 11 # how many labels to see on axis x
        step_x = int(nx / (no_labels - 1)) # step between consecutive labels
        x_positions = np.arange(3,nx+1,step_x) # pixel count at label position
        x_labels = x[3::step_x] # labels you want to see
        plt.xticks(x_positions, x_labels)
        plt.yticks(x_positions, x_labels)
        plt.title('b)', fontsize=20, fontweight="bold") # Edit title
        plt.xlabel('Participant No.', fontsize=18)
        plt.ylabel('Participant No.', fontsize=18)
        plt.savefig('/dhcp/fmri_anna_graham/GKgit/head_position/FINGER/finger_figures/predFC_figures/RSM_predFC_older.jpg') #edit older/younger here!
    if youngersession:
        plt.figure(figsize=(12,8))
        c = plt.imshow(comparesess, cmap='Blues') # Edit color here
        plt.colorbar(c)
        x = np.arange(1,45,1) # the grid to which your data corresponds
        nx = x.shape[0]
        no_labels = 11 # how many labels to see on axis x
        step_x = int(nx / (no_labels - 1)) # step between consecutive labels
        x_positions = np.arange(3,nx+1,step_x) # pixel count at label position
        x_labels = x[3::step_x] # labels you want to see
        plt.xticks(x_positions, x_labels)
        plt.yticks(x_positions, x_labels)
        plt.title('a)', fontsize=20, fontweight="bold") # Edit title
        plt.xlabel('Participant No.', fontsize=18)
        plt.ylabel('Participant No.', fontsize=18)
        plt.savefig('/dhcp/fmri_anna_graham/GKgit/head_position/FINGER/finger_figures/predFC_figures/RSM_predFC_younger.jpg') #edit older/younger here!

    # Identifying whether within_predFC is the highest correlation:
    match = np.diag(comparesess) == np.max(comparesess, axis=0)
    for ind in range(44):
        print('%s'%(match[ind]))
    print()   


    ## Calculating the difference in predFC correlation meanwithin vs meanbetween:
    iu1 = np.triu_indices(nsubj, k=1)
    within = np.diag(comparesess)
    between = comparesess[iu1]

    mean_within = np.mean(within)
    print(f'Second level correlation - Mean within subject:{mean_within}')
    mean_between = np.mean(between)
    print(f'Second level correlation - Mean between subject:{mean_between}')
    mean_within_minus_between = mean_within - mean_between
    print(f'Second Level Correlation - Mean within subject minus between subject:{mean_within_minus_between}')
    print()



if __name__ == '__main__':
    
    ##SWITCHES for choosing acrosssessions or younger_vs_older session:
    oldersession = False
    youngersession = False
    acrosssession = True

    allsnr_iu1, nsubj, nsess, nhp, nedges = loadarray()
    predicted = predicted_fc(allsnr_iu1, nsubj, nsess, nhp, nedges)
    fingerprint(oldersession, youngersession, acrosssession, predicted, nsubj)

