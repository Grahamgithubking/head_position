import pandas as pd
import numpy as np
import glob
import os
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
from scipy import stats
import seaborn as sns


def loadfc():
    ## Choose either of these:
    # allfc = np.load('/dhcp/fmri_anna_graham/GKgit/finger_npy/fc_416.npy')
    allfc = np.load ('/dhcp/fmri_anna_graham/GKgit/finger_npy/fc_96.npy')
    
    
    sz=allfc.shape
    print(f"starting allfc shape is: {sz}")
    nsubj=sz[0]
    nsess=sz[1]
    nroi=sz[2]
    allfc_reshaped=np.reshape(allfc, (nsubj*nsess, nroi, nroi))
    print('The shape of allfc_reshaped is:')
    print(allfc_reshaped.shape)

    return allfc_reshaped, nsubj, nsess, nroi

def loadsnr(snrproduct, nsubj, nsess, nroi):
    ## Choose either of these:
    # allsnr = np.load('/dhcp/fmri_anna_graham/GKgit/snr_npy/416_snr_416_erode1.npy')
    # allsnr = np.load('/dhcp/fmri_anna_graham/GKgit/snr_npy/416_snr_true.npy')
    allsnr = np.load('/dhcp/fmri_anna_graham/GKgit/snr_npy/96_snr_true.npy')
    

    mz=allsnr.shape
    print(f"The shape of allsnr is: {mz}")

    allsnr_reshaped=np.reshape(allsnr, (nsubj*nsess, nroi))
    print('The shape of allsnr_reshaped is:')
    print(allsnr_reshaped.shape)
    print()

    print('allsnr_reshaped at [:3, 153] are:')
    print(allsnr_reshaped[:3,153])
    print()

    ### Multiplying the SNR across all edges for each ROI:
    if snrproduct:
        snrbysnr=np.zeros((nsubj*nsess,nroi,nroi))
        for ind in range(nsubj*nsess):
            snr=allsnr_reshaped[ind,:]
            snrbysnr[ind,:,:]=np.outer(snr, snr) #instead of using snr*snr.T
        print('The shape of snrbysnr is:')
        print(snrbysnr.shape)
        print()
        allsnr_reshaped=snrbysnr

    
    return allsnr_reshaped


def select_product(item_reshaped, regions, nsubj, nsess):
    # To use with 3-d arrays (where the 2nd and 3rd dimensions are the edge correlation/product)
    mat = np.empty([nsubj*nsess,len(regions),len(regions)])
    for subind, subj in enumerate(np.arange(nsubj*nsess)):
        for inda, regiona in enumerate(regions):
            for indb, regionb in enumerate(regions):
                product_temp = item_reshaped[subj,regiona,regionb]
                mat[subind,inda,indb] = product_temp
    print(f"These are some values of select_product: {mat[:3,:,:]}")
    print()
    iu1=np.triu_indices(len(regions), k=1) #selecting the upper triangle of a matrix (not including main diagonal)
    mat_iu1=np.array([x[iu1]for x in mat])
    print(f"The shape of mat_iu1 is: {mat_iu1.shape}")
    print(f"These are some items of mat_iu1: {mat_iu1[:3,:]}")
    print()

    return mat_iu1

def selectsnr_adding(item_reshaped, regions, nsubj, nsess):
    # To use with 2-d arrays (where roi values are added to get an edge value)
    numedges = int(((len(regions))-1)*((len(regions))/2))
    mat = np.empty([nsubj*nsess, numedges])
    for subind, subj in enumerate(np.arange(nsubj*nsess)):
        for inda, regiona in enumerate(regions):
            for indb, regionb in enumerate(regions):
                if (indb > inda): # Choosing unique edges and avoiding any previously added snr pairs
                    addsnr = item_reshaped[subj,regiona] + item_reshaped[subj,regionb]
                    print(addsnr)

                    ###### ???? how to add addsnr value into numpy array 'mat'[nsubj*nsess, numedges] ??????
                        
    print(f"These are some values of selectsnr_adding: {mat[:3,:]}")
    print()

    return mat


if __name__ == '__main__':

    #SWITCH:
    snrproduct = True  # A switch to decide if using [product vs addition] of SNR values to get the edge

    allfc_reshaped, nsubj, nsess, nroi = loadfc()
    allsnr_reshaped = loadsnr(snrproduct,nsubj,nsess,nroi)

    # 154	17Networks_LH_DefaultA_pCunPCC_1
    # 161	17Networks_LH_DefaultA_PFCm_1
    # 167	17Networks_LH_DefaultB_Temp_1
    # 175	17Networks_LH_DefaultB_PFCd_1
    # 363	17Networks_RH_DefaultA_pCunPCC_1
    # 368	17Networks_RH_DefaultA_PFCm_1
    # 374	17Networks_RH_DefaultB_Temp_1
    # 377	17Networks_RH_DefaultB_PFCd_1

    roisonebased = list(range(149,195)) # all DMN on LH Side (range up to 194)
    # roisonebased = [154,161,167,175,363,368,374,377] # Bilateral LH and RH
    # roisonebased = [154,161,167,175] # 4roi of DMN on LH Side

    roiszerobased=[x-1 for x in roisonebased]
    print(f"These are the roiszerobased: {roiszerobased}")
    print(f"The length of roiszerobased is: {len(roiszerobased)}")
    print()

    ### Finding out the ROI pair with highest meanFC across sessions:
    roiconn = []
    for roiA in roiszerobased:
        for roiB in roiszerobased:
            if not roiA==roiB:
                fc_all = allfc_reshaped[:,roiA,roiB]
                fc_mean = np.mean(fc_all)
                roiconn.append([roiA,roiB,fc_mean])
    df = pd.DataFrame(roiconn, columns =['RoiA', 'RoiB', 'meanFC'])
    print(df)
    print()
    df_sorted = df.sort_values(by = 'meanFC')
    print(f"This is df_sorted by meanFC:")
    print(df_sorted)






    ### Enter the selected "zerobased" pair in here:
    ### 17Networks_LH_DefaultB_IPL_2	205	62	85	0
    ### 17Networks_LH_DefaultB_PFCd_1	205	63	77	0
    pairzerobased = [173,174]

    ### Calculating across all sessions:
    fc_iu1 = select_product(allfc_reshaped, pairzerobased, nsubj, nsess)
    if snrproduct:
        ### If using select_product:
        snr_iu1 = select_product(allsnr_reshaped, pairzerobased, nsubj, nsess)
    else:
        ### If using selectsnr_adding:
        snr_iu1 = selectsnr_adding(allsnr_reshaped, pairzerobased, nsubj, nsess)

    # #######  
    # Option: Standardization?? I think not!!!!!!
    # fc_iu1 = stats.zscore(fc_iu1, axis=1)
    # snr_iu1 = stats.zscore(snr_iu1, axis=1)

    
    ### Plotting:
    plt.figure(1)
    plt.scatter(np.ravel(fc_iu1), np.ravel(snr_iu1))
    plt.title('96 sessions \n LH_DefaultB_IPL_2 \n LH_DefaultB_PFCd_1')
    plt.ylabel('SNR-individual')
    plt.xlabel('FC')
    plt.savefig('/dhcp/fmri_anna_graham/GKgit/head_position/FINGER/finger_figures/fingerfigures_misc/SNR_covary_FC.jpg')

    plt.figure(2)
    sns.kdeplot(x=np.ravel(fc_iu1), y=np.ravel(snr_iu1))
    plt.title('96 sessions \n LH_DefaultB_IPL_2 \n LH_DefaultB_PFCd_1')
    plt.ylabel('SNR-individual')
    plt.xlabel('FC')
    plt.savefig('/dhcp/fmri_anna_graham/GKgit/head_position/FINGER/finger_figures/fingerfigures_misc/SNR_covary_FC_kde.jpg')

