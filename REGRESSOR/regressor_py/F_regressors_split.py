## To be continued....

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


def loadsnr(snrtrue):
    # Loading SNR values for two splits per session
    # for all sessions of the 48 particpants with both preterm and term sessions
    if snrtrue:
        snr_all = np.load('/dhcp/fmri_anna_graham/GKgit/snr_npy/192_snrtrue.npy')
    else:
        snr_all = np.load('/dhcp/fmri_anna_graham/GKgit/snr_npy/192hp_snr.npy')
    mz=snr_all.shape
    print('starting snr_all shape is:')
    print(mz) #snr_all has shape [48, 2, 2, 387]
    print()
    nsubj=mz[0]
    nsess=mz[1]
    nhp=mz[2]
    nroi=mz[3]
    print('Some values of snr_all are:')
    print(snr_all[:2,:1,:1,:10])
    print()

    snr_all_reshaped=np.reshape(snr_all, (nsubj*nsess*nhp, nroi))
    print('The shape of snr_all_reshaped is:')
    print(snr_all_reshaped.shape) # snr_all_reshaped has shape [192, 387]
    print()

    ## Products/edges of SNR values (to match the format of the original snr products used to deduce the OLS parameters c and m):
    snrbysnr=np.zeros((nsubj*nsess*nhp, nroi, nroi))
    for ind in range(nsubj*nsess*nhp):
        snr=snr_all_reshaped[ind,:]
        snrbysnr[ind,:,:]=np.outer(snr, snr) #instead of using snr*snr.T
    print('The shape of snrbysnr is:')
    print(snrbysnr.shape) # This is [192,387,387]
    print()
    iu1 = np.triu_indices(nroi, k=1)
    allsnr_iu1=np.array([x[iu1]for x in snrbysnr]) #for each session, selecting the upper triangle of the 387x387 matrix
    print('The shape of allsnr_iu1 is:')
    print(allsnr_iu1.shape)
    print()
    print('Some values of allsnr_iu1 are:')
    print(allsnr_iu1[:2,:5])
    print()
    nedges=allsnr_iu1.shape[1]
    print('The number of edges is:')
    print(nedges)
    print()

    allsnr_sz=allsnr_iu1
    # ### Standardizing the SNR edge values:
    # smean=np.mean(allsnr_iu1, axis=1)
    # sstd=np.std(allsnr_iu1, axis=1)
    # allsnr_sz=np.zeros((nsubj*nsess*nhp,nedges))
    # for ind in range(nsubj*nsess*nhp):
    #     allsnr_sz[ind,:]=(allsnr_iu1[ind,:]-smean[ind,])/sstd[ind,]
    # print('Some values of allsnr_sz are:')
    # print(allsnr_sz[:2,25200:25210])
    # print()

    return allsnr_sz, nedges

def loadfc():
    allfc = np.load('/dhcp/fmri_anna_graham/GKgit/finger_npy/192fc.npy')
    sz=allfc.shape
    print('Starting allfc shape is:')
    print(sz)
    nsubj=sz[0]
    nsess=sz[1]
    nsplit=sz[2]
    nroi=sz[3]

    ## Reshaping 
    allfc_reshaped=np.reshape(allfc, (nsubj*nsess*nsplit, nroi, nroi))
    print('The shape of allfc_reshaped is:')
    print(allfc_reshaped.shape)

    iu1 = np.triu_indices(nroi, k=1)
    allfc_iu1=np.array([x[iu1]for x in allfc_reshaped]) #for each session, selecting the upper triangle of the 387x387 matrix
    print('The shape of allfc_iu1 is:')
    print(allfc_iu1.shape)
    print()
    print('Some values of allfc_iu1 are:')
    print(allfc_iu1[:2,:5])
    print()

    # ### Standardising the FC values:
    # allfc_iu1=stats.zscore(allfc_iu1, axis=1)

    return allfc_iu1, nsubj


def regressors(snrtrue,sess1,sess2,splitA,splitB,snrfcya,snrfcyb,nsubj,nedges,allsnr_sz,allfc_iu1,orthog,subjectlevel,edgelevel):

    if sess1:
        if splitA:
            # Selecting SNR_Ss1_SpA:
            snr=allsnr_sz[0::4,:]
            # Selecting FC_Ss1_SpA:
            fcx=allfc_iu1[0::4,:]
            # Selecting FC_Ss1_SpB:
            fcy=allfc_iu1[1::4,:]
        if snrfcya:
            snr=allsnr_sz[0::4,:]
            fcx=allfc_iu1[1::4,:]
            fcy=allfc_iu1[0::4,:]
        if splitB:
            # Selecting SNR_Ss1_SpB:
            snr=allsnr_sz[1::4,:]
            # Selecting FC_Ss1_SpB:
            fcx=allfc_iu1[1::4,:]
            # Selecting FC_Ss1_SpA:
            fcy=allfc_iu1[0::4,:]
        if snrfcyb:
            snr=allsnr_sz[1::4,:]
            fcx=allfc_iu1[0::4,:]
            fcy=allfc_iu1[1::4,:]
        
    if sess2:
        if splitA:
            # Selecting SNR_Ss2_SpA:
            snr=allsnr_sz[2::4,:]
            # Selecting FC_Ss2_SpA:
            fcx=allfc_iu1[2::4,:]
            # Selecting FC_Ss2_SpB:
            fcy=allfc_iu1[3::4,:]
        if snrfcya:
            snr=allsnr_sz[2::4,:]
            fcx=allfc_iu1[3::4,:]
            fcy=allfc_iu1[2::4,:]
        if splitB:
            # Selecting SNR_Ss2_SpB:
            snr=allsnr_sz[3::4,:]
            # Selecting FC_Ss2_SpB:
            fcx=allfc_iu1[3::4,:]
            # Selecting FC_Ss2_SpA:
            fcy=allfc_iu1[2::4,:]
        if snrfcyb:
            snr=allsnr_sz[3::4,:]
            fcx=allfc_iu1[2::4,:]
            fcy=allfc_iu1[3::4,:]

    if edgelevel:

        parameters=np.zeros((nedges,3)) # <3> parameters here as yes/no constants are selected!
        p=np.zeros(nedges)
        r2=np.zeros(nedges)

        ### Standardising the values:
        snr=stats.zscore(snr, axis=1)
        fcx=stats.zscore(fcx, axis=1)
        fcy=stats.zscore(fcy, axis=1)

        for r in range(nedges):
            if orthog:
                # Taking out part of fcx that looks like snr, aka orthogonalisation:
                fcx[:,r] = fcx[:,r]-(np.dot(fcx[:,r], snr[:,r])/np.dot(snr[:,r], snr[:,r]))*snr[:,r]
            
            # Multiple linear regression - independent varibles SNR and FCx
            x=np.stack((snr[:,r],fcx[:,r]), axis=0)
            y=fcy[:,r]

            X2 = sm.add_constant(x.T)
            results = sm.OLS(y, X2).fit()
            ## y = c + mx + bz where (c, m, b) are parameters
            ## NOTE: param1 is the constant_value, param2 is the 1st coefficient, param2 is the 2nd coefficient
            parameters[r] = results.params
            p[r] = results.pvalues[1]
            r2[r] = results.rsquared
        print('The parameters array size is:')
        print(parameters.shape)
        print()
        print('Some parameters are:')
        print(parameters[:10,:])
        print()

        if snrtrue:
            if sess1:
                if splitA:
                    if orthog:
                        np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/74691_true_orthog1A.npy', parameters)
                    else:
                        np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/74691_true_1A.npy', parameters)
                if splitB:
                    if orthog:
                        np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/74691_true_orthog1B.npy', parameters)
                    else:
                        np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/74691_true_1B.npy', parameters)
            if sess2:
                if splitA:
                    if orthog:
                        np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/74691_true_orthog2A.npy', parameters)
                    else:
                        np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/74691_true_2A.npy', parameters)
                if splitB:
                    if orthog:
                        np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/74691_true_orthog2B.npy', parameters)
                    else:
                        np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/74691_true_2B.npy', parameters)
        else:
            if sess1:
                if splitA:
                    if orthog:
                        np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/74691_params_orthog1A.npy', parameters)
                    else:
                        np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/74691_params_1A.npy', parameters)
                if splitB:
                    if orthog:
                        np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/74691_params_orthog1B.npy', parameters)
                    else:
                        np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/74691_params_1B.npy', parameters)
                if snrfcya:
                    if orthog:
                        np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/74691_params_orthog1snrfcya.npy', parameters)
                    else:
                        np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/74691_params_1snrfcya.npy', parameters)
                if snrfcyb:
                    if orthog:
                        np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/74691_params_orthog1snrfcyb.npy', parameters)
                    else:
                        np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/74691_params_1snrfcyb.npy', parameters)
            if sess2:
                if splitA:
                    if orthog:
                        np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/74691_params_orthog2A.npy', parameters)
                    else:
                        np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/74691_params_2A.npy', parameters)
                if splitB:
                    if orthog:
                        np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/74691_params_orthog2B.npy', parameters)
                    else:
                        np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/74691_params_2B.npy', parameters)
                if snrfcya:
                    if orthog:
                        np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/74691_params_orthog2snrfcya.npy', parameters)
                    else:
                        np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/74691_params_2snrfcya.npy', parameters)
                if snrfcyb:
                    if orthog:
                        np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/74691_params_orthog2snrfcyb.npy', parameters)
                    else:
                        np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/74691_params_2snrfcyb.npy', parameters)


    if subjectlevel:

        parameters=np.zeros((nsubj,3)) # <3> parameters here as yes/no constants are selected!
        p=np.zeros(nsubj)
        r2=np.zeros(nsubj)

        snr=stats.zscore(snr, axis=1)
        fcx=stats.zscore(fcx, axis=1)
        fcy=stats.zscore(fcy, axis=1)

        for r in range(nsubj):
            if orthog:
                # Taking out part of fcx that looks like snr, aka orthogonalisation:
                fcx[r,:] = fcx[r,:]-(np.dot(fcx[r,:], snr[r,:])/np.dot(snr[r,:], snr[r,:]))*snr[r,:]
            
            # Multiple linear regression - independent varibles SNR and FCx
            x=np.stack((snr[r,:],fcx[r,:]), axis=0)
            y=fcy[r,:]

            X2 = sm.add_constant(x.T)
            results = sm.OLS(y, X2).fit()
            ## y = c + mx + bz where (c, m, b) are parameters
            ## NOTE: param1 is the constant_value, param2 is the 1st coefficient, param2 is the 2nd coefficient
            parameters[r] = results.params
            p[r] = results.pvalues[1]
            r2[r] = results.rsquared
        print('The parameters array size is:')
        print(parameters.shape)
        print()
        print('Some parameters are:')
        print(parameters[:12,:])
        print()

        if snrtrue:
            if sess1:
                if splitA:
                    if orthog:
                        np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/48_true_orthog1A.npy', parameters)
                    else:
                        np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/48_true_1A.npy', parameters)
                if splitB:
                    if orthog:
                        np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/48_true_orthog1B.npy', parameters)
                    else:
                        np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/48_true_1B.npy', parameters)
            if sess2:
                if splitA:
                    if orthog:
                        np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/48_true_orthog2A.npy', parameters)
                    else:
                        np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/48_true_2A.npy', parameters)
                if splitB:
                    if orthog:
                        np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/48_true_orthog2B.npy', parameters)
                    else:
                        np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/48_true_2B.npy', parameters)
        else:
            if sess1:
                if splitA:
                    if orthog:
                        np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/48_params_orthog1A.npy', parameters)
                    else:
                        np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/48_params_1A.npy', parameters)
                if splitB:
                    if orthog:
                        np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/48_params_orthog1B.npy', parameters)
                    else:
                        np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/48_params_1B.npy', parameters)
            if sess2:
                if splitA:
                    if orthog:
                        np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/48_params_orthog2A.npy', parameters)
                    else:
                        np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/48_params_2A.npy', parameters)
                if splitB:
                    if orthog:
                        np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/48_params_orthog2B.npy', parameters)
                    else:
                        np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/48_params_2B.npy', parameters)

if __name__ == '__main__':

    #SWITCHES:
    
    subjectlevel=False
    edgelevel=True

    orthog=True  #  Turn this of off if using snrfcy!

    snrtrue=False  # If using true snr (True) versus Snr_Coil (False)

    sess1=False
    sess2=True

    splitA=True  # Where SNR (Coil/True) is from a 'different' split as FCy
    splitB=False

    snrfcya=False  # Where SNR-Coil is from the 'SAME' split as FCy
    snrfcyb=False


    allsnr_sz, nedges = loadsnr(snrtrue)
    allfc_iu1, nsubj = loadfc()
    regressors(snrtrue,sess1,sess2,splitA,splitB,snrfcya,snrfcyb,nsubj,nedges,allsnr_sz,allfc_iu1,orthog,subjectlevel,edgelevel)

