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


def loadsnr(snrcoil, snrtrue):
   
    def reshapesnr(snr_all):
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

        return allsnr_iu1

    if snrcoil:
        snrcoil_all = np.load('/dhcp/fmri_anna_graham/GKgit/snr_npy/192hp_snr.npy')
        allsnrcoil_iu1=reshapesnr(snrcoil_all)
    if snrtrue:
        snrtrue_all = np.load('/dhcp/fmri_anna_graham/GKgit/snr_npy/192hp_snrtrue.npy')
        allsnrtrue_iu1=reshapesnr(snrtrue_all)

    return allsnrcoil_iu1, allsnrtrue_iu1

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

    nedges=allfc_iu1.shape[1]
    print('The number of edges is:')
    print(nedges)
    print()

    return allfc_iu1, nsubj, nedges


def regressors(subjectlevel,edgelevel,orthog,snrcoil,snrtrue,sess1,sess2,splitA,splitB,allsnrcoil_iu1,allsnrtrue_iu1,allfc_iu1, nsubj, nedges):

    def defineregress(allsnr_iu1,allfc_iu1):
        if sess1:
            if splitA:
                # Selecting SNR_Ss1_SpA:
                snr=allsnr_iu1[0::4,:]
                # Selecting FC_Ss1_SpA:
                fcx=allfc_iu1[0::4,:]
                # Selecting FC_Ss1_SpB:
                fcy=allfc_iu1[1::4,:]
            if splitB:
                # Selecting SNR_Ss1_SpB:
                snr=allsnr_iu1[1::4,:]
                # Selecting FC_Ss1_SpB:
                fcx=allfc_iu1[1::4,:]
                # Selecting FC_Ss1_SpA:
                fcy=allfc_iu1[0::4,:]
        if sess2:
            if splitA:
                # Selecting SNR_Ss2_SpA:
                snr=allsnr_iu1[2::4,:]
                # Selecting FC_Ss2_SpA:
                fcx=allfc_iu1[2::4,:]
                # Selecting FC_Ss2_SpB:
                fcy=allfc_iu1[3::4,:]
            if splitB:
                # Selecting SNR_Ss2_SpB:
                snr=allsnr_iu1[3::4,:]
                # Selecting FC_Ss2_SpB:
                fcx=allfc_iu1[3::4,:]
                # Selecting FC_Ss2_SpA:
                fcy=allfc_iu1[2::4,:]

        return snr, fcx, fcy
    
    coilsnr, fcx, fcy = defineregress(allsnrcoil_iu1, allfc_iu1)
    truesnr, fcx, fcy = defineregress(allsnrtrue_iu1, allfc_iu1)

    if edgelevel:

        parameters=np.zeros((nedges,4)) # <4> parameters here as yes/no constants are selected!

        ### Standardising the values:
        coilsnr=stats.zscore(coilsnr, axis=1)
        truesnr=stats.zscore(truesnr, axis=1)
        fcx=stats.zscore(fcx, axis=1)
        fcy=stats.zscore(fcy, axis=1)

        for r in range(nedges):
            if orthog:
                # Taking out part of fcx that looks like truesnr, aka orthogonalisation:
                fcx[:,r] = fcx[:,r]-(np.dot(fcx[:,r], truesnr[:,r])/np.dot(truesnr[:,r], truesnr[:,r]))*truesnr[:,r]
        
            x=np.stack((coilsnr[:,r],truesnr[:,r],fcx[:,r]), axis=0)
            y=fcy[:,r]

            X2 = sm.add_constant(x.T)
            results = sm.OLS(y, X2).fit()
            ## y = c + mx + bz where (c, m, b) are parameters
            ## NOTE: param1 is the constant_value, param2 is the 1st coefficient, param2 is the 2nd coefficient
            parameters[r] = results.params
        
        print('The parameters array size is:')
        print(parameters.shape)
        print()
        print('Some parameters are:')
        print(parameters[:10,:])
        print()

    
        if sess1:
            if splitA:
                if orthog:
                    np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/74691_snrboth_orthog1A.npy', parameters)
                else:
                    np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/74691_snrboth_1A.npy', parameters)
            if splitB:
                if orthog:
                    np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/74691_snrboth_orthog1B.npy', parameters)
                else:
                    np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/74691_snrboth_1B.npy', parameters)
        if sess2:
            if splitA:
                if orthog:
                    np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/74691_snrboth_orthog2A.npy', parameters)
                else:
                    np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/74691_snrboth_2A.npy', parameters)
            if splitB:
                if orthog:
                    np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/74691_snrboth_orthog2B.npy', parameters)
                else:
                    np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/74691_snrboth_2B.npy', parameters)
        

if __name__ == '__main__':

    #SWITCHES:
    
    subjectlevel=False
    edgelevel=True

    orthog=True

    snrcoil=True
    snrtrue=True 

    sess1=False
    sess2=True

    splitA=True  # Where SNR (Coil/True) is from a 'different' split as FCy
    splitB=False

    allsnrcoil_iu1, allsnrtrue_iu1 = loadsnr(snrcoil, snrtrue)
    allfc_iu1, nsubj, nedges = loadfc()
    regressors(subjectlevel,edgelevel,orthog,snrcoil,snrtrue,sess1,sess2,splitA,splitB,allsnrcoil_iu1,allsnrtrue_iu1,allfc_iu1, nsubj, nedges)

