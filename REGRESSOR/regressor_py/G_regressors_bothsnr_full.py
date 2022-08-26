#### All original code is copyright of Graham King, Trinity College Institute of Neuroscience.
### Using predicted SNR-HeadPosition and corresponding Fuctional Connectivity data to calculate the parameters of 74,691 linear regressors across 416 participants with one session.
### Note: The raw SNR data from these 416 participants was used to create the SNR-Coil (SNR of the Coil)

import pandas as pd
import numpy as np
import glob
import os
import matplotlib.pyplot as plt
from scipy import stats
import pingouin as pg
import statsmodels.api as sm
import seaborn as sns


def loadfc(one_sessid, two_sessid):
    #Loading and Reshaping both predicted SNR-HeadPosition data and FC data
    
    if one_sessid:
        allfc=np.load('/dhcp/fmri_anna_graham/GKgit/finger_npy/fc_416.npy')
        sz=allfc.shape
        print('one_sessid allfc shape is:')
        print(sz)
        nsubj=sz[0]
        nsess=sz[1]
        nroi=sz[2]
    if two_sessid:
        allfc=np.load('/dhcp/fmri_anna_graham/GKgit/finger_npy/fc_96.npy')
        sz=allfc.shape
        print('two_sessid allfc shape is:')
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

    ###### Getting the Fishers r2z standardization:
    # allfc_fish=np.arctanh(allfc_iu1)


    return allfc_iu1, nsubj, nsess, nroi, nedges

def loadsnr(one_sessid, two_sessid, snrcoil, snrtrue, nsubj, nsess, nroi, nedges):

    def reshapesnr(snr_all):
        print('The shape of snr_all is:')
        print(snr_all.shape)
        print()
        snr_all_reshaped=np.reshape(snr_all, (nsubj*nsess, nroi))
        print('The shape of snr_all_reshaped is:')
        print(snr_all_reshaped.shape)
        print()

        # Multiplying the SNR across all edges for each ROI:
        snrbysnr=np.zeros((nsubj*nsess,nroi,nroi))
        for ind in range(nsubj*nsess):
            snr=snr_all_reshaped[ind,:]
            snrbysnr[ind,:,:]=np.outer(snr, snr) #instead of using snr*snr.T
        print('The shape of snrbysnr is:')
        print(snrbysnr.shape)
        print()

        iu1=np.triu_indices(nroi, k=1) #selecting the upper triangle of a matrix
        allsnr_iu1=np.array([x[iu1]for x in snrbysnr]) #for each session, selecting the upper triangle of the 387x387 matrix
        print('The shape of allsnr_iu1 is:')
        print(allsnr_iu1.shape)
        print('Some values of allsnr_iu1 are:')
        print(allsnr_iu1[:2,:5])
        print()

        return allsnr_iu1

    if one_sessid:
        snrone_all=np.load('/dhcp/fmri_anna_graham/GKgit/snr_npy/416_snr_416_erode1.npy')
        allsnrone_iu1=reshapesnr(snrone_all)
    if two_sessid:
        if snrcoil:
            snrcoil_all=np.load('/dhcp/fmri_anna_graham/GKgit/snr_npy/96_snr_416_erode1.npy')
            allsnrcoil_iu1=reshapesnr(snrcoil_all)
        if snrtrue:
            snrtrue_all=np.load('/dhcp/fmri_anna_graham/GKgit/snr_npy/96_snr_true.npy')
            allsnrtrue_iu1=reshapesnr(snrtrue_all)
         
    return allsnrcoil_iu1, allsnrtrue_iu1


def regressors(one_sessid, two_sessid, allfc_iu1, allsnrcoil_iu1, allsnrtrue_iu1, nedges, nsubj, snrcoil, snrtrue, sess1, sess2, subjectlevel, edgelevel, orthog):

    if one_sessid:
        parameters=np.zeros((nedges,2)) #nedges OR specify number nedge of 74,691 here!
        p=np.zeros(nedges)
        r2=np.zeros(nedges)

        for r in range(nedges):
            x=allsnrone_iu1[:,r]
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
        np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/74691_params_snr_416_erode1.npy', parameters)
        np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/74691_p_snr_416_erode1.npy', p)
        np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/74691_r2_snr_416_erode1.npy', r2)

    if two_sessid: # allsnr_iu1 has size [96,74691]  # allfc_iu1 has size [96,74691]
        
        def defineregress(allsnr_iu1,allfc_iu1):
            if sess1:
                snr=allsnr_iu1[0::2,:]
                fcx=allfc_iu1[0::2,:]
                fcy=allfc_iu1[1::2,:]
            if sess2:
                snr=allsnr_iu1[1::2,:]
                fcx=allfc_iu1[1::2,:]
                fcy=allfc_iu1[0::2,:]

            return snr, fcx, fcy
        
        coilsnr, fcx, fcy =defineregress(allsnrcoil_iu1, allfc_iu1)
        truesnr, fcx, fcy =defineregress(allsnrtrue_iu1, allfc_iu1)


        if edgelevel:
            parameters=np.zeros((nedges,4)) # <4> parameters here as <yes> constants are selected!

            ### Standardising the SNR and FC values:
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
            # DF = pd.DataFrame(parameters)
            # DF.to_csv("/dhcp/fmri_anna_graham/GKgit/head_position/REGRESSOR/regressor_py/parameters_csv.csv")

            if sess1:
                if orthog:
                    np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/74691_params_1_snrboth_orthog.npy', parameters)
                else:
                    np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/74691_params_1_snrboth.npy', parameters)
            if sess2:
                if orthog:
                    np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/74691_params_2_snrboth_orthog.npy', parameters)
                else:
                    np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/74691_params_2_snrboth.npy', parameters)

    
if __name__ == '__main__':

    one_sessid=False
    two_sessid=True

    subjectlevel=False  #subjectlevel is: regressors across edges
    edgelevel=True  #edgelevel is: regressors across subjects

    orthog=True # Turn this off if using snrfcy!

    snrcoil=True  # aka SNR-Coil
    snrtrue=True # Turn this off if using snrfcy!

    sess1=False
    sess2=True
    
    allfc_iu1, nsubj, nsess, nroi, nedges = loadfc(one_sessid, two_sessid)
    allsnrcoil_iu1, allsnrtrue_iu1 = loadsnr(one_sessid, two_sessid, snrcoil, snrtrue, nsubj, nsess, nroi, nedges)
    regressors(one_sessid, two_sessid, allfc_iu1, allsnrcoil_iu1, allsnrtrue_iu1, nedges, nsubj, snrcoil, snrtrue, sess1, sess2, subjectlevel, edgelevel, orthog)


 
