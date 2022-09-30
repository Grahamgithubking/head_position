import pandas as pd
import numpy as np
import glob
import os
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
from scipy import stats
import seaborn as sns
import pingouin as pg


def loadfc(allfc):
    sz=allfc.shape
    print(f"starting allfc shape is: {sz}")
    nsubj=sz[0]
    nsess=sz[1]
    nroi=sz[2]
    allfc_reshaped=np.reshape(allfc, (nsubj*nsess, nroi, nroi))
    print('The shape of allfc_reshaped is:')
    print(allfc_reshaped.shape)

    return allfc_reshaped, nsubj, nsess, nroi

def loadsnr(allsnr, snrproduct, nsubj, nsess, nroi):
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


if __name__ == '__main__':

    #SWITCH:
    snrproduct = True  # A switch to decide if using [product vs addition] of SNR values to get the edge

    ### Choose either of these:
    # allfc416 = np.load('/dhcp/fmri_anna_graham/GKgit/finger_npy/fc_416.npy')
    allfc96 = np.load ('/dhcp/fmri_anna_graham/GKgit/finger_npy/fc_96.npy')


    ### Choose which of these:
    # allsnr = np.load('/dhcp/fmri_anna_graham/GKgit/snr_npy/416_snr_416_erode1.npy')
    # allsnr = np.load('/dhcp/fmri_anna_graham/GKgit/snr_npy/416_snr_true.npy')
    allsnrindividual = np.load('/dhcp/fmri_anna_graham/GKgit/snr_npy/96_snr_true.npy')
    allsnrgroup = np.load('/dhcp/fmri_anna_graham/GKgit/snr_npy/96_snr_416_erode1.npy')


    allfc_reshaped, nsubj, nsess, nroi = loadfc(allfc96)

    allsnr_reshaped_indiv = loadsnr(allsnrindividual,snrproduct,nsubj,nsess,nroi)
    allsnr_reshaped_group = loadsnr(allsnrgroup,snrproduct,nsubj,nsess,nroi)

    roisonebased = list(range(142,188)) # all 46roi DMN on LH Side, allowing for fact that 7 rois were trimmed prior to 400roi 149-194+1!!!!!!!
    # roisonebased = list(range(345,378)) # all 33roi DMN on RH Side, allowing for fact that 13 rois were trimmed prior to 400roi 358-390+1!!!!

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


    ### Enter the selected "ONEbased" pair in here:
    # 174	17Networks_LH_DefaultB_IPL_2	205	62	85	0
    # 175	17Networks_LH_DefaultB_PFCd_1	205	63	77	0
    # 372	17Networks_RH_DefaultA_PFCm_5	249	255	5	0
    # 373	17Networks_RH_DefaultA_PFCm_6	249	255	6	0

    pairzerobased = [173,174] # LH Side, fc 0.5749
    # pairzerobased = [371,372] # RH Side, fc 0.5642


    ### Calculating across all sessions:
    fc_iu1 = select_product(allfc_reshaped, pairzerobased, nsubj, nsess)
    if snrproduct:
        ### If using select_product:
        snr_iu1_indiv = select_product(allsnr_reshaped_indiv, pairzerobased, nsubj, nsess)
        snr_iu1_group = select_product(allsnr_reshaped_group, pairzerobased, nsubj, nsess)
    
    print(fc_iu1.shape)
    print(snr_iu1_group.shape)

    # #######  
    # Option: Standardization?? I think not!!!!!!
    # fc_iu1 = stats.zscore(fc_iu1, axis=1)
    # snr_iu1 = stats.zscore(snr_iu1, axis=1)

    
    ### Plotting:
    plt.figure(1)
    plt.scatter(x=(fc_iu1), y=(snr_iu1_indiv), c='red')
    # plt.title('?? sessions \n TBC')
    plt.ylabel('SNR-individual')
    plt.xlabel('FC')
    plt.legend(['r=0.248, p<0.015'], loc='lower right', fontsize=13)
    plt.savefig('/dhcp/fmri_anna_graham/GKgit/head_position/FINGER/finger_figures/fingerfigures_misc/SNRcovaryFC_snrindiv.jpg')

    pc1 = pg.corr(np.ravel(fc_iu1), np.ravel(snr_iu1_indiv), tail='two-sided', method='pearson')
    print('The pearson correlation coefficient for FC and SNR-individual is:')
    print(pc1)
    print()

    plt.figure(2)
    plt.scatter(x=(fc_iu1), y=(snr_iu1_group), c='blue')
    # plt.title('?? sessions \n TBC')
    plt.ylabel('SNR-group')
    plt.xlabel('FC')
    plt.legend(['r=0.253, p<0.013'], loc='lower right', fontsize=13)
    plt.savefig('/dhcp/fmri_anna_graham/GKgit/head_position/FINGER/finger_figures/fingerfigures_misc/SNRcovaryFC_snrgroup.jpg')

    pc2 = pg.corr(np.ravel(fc_iu1), np.ravel(snr_iu1_group), tail='two-sided', method='pearson')
    print('The pearson correlation coefficient for FC and SNR-group is:')
    print(pc2)
    print()


    plt.figure(3)
    sns.kdeplot(x=np.ravel(fc_iu1), y=np.ravel(snr_iu1_indiv), color='red')
    sns.kdeplot(x=np.ravel(fc_iu1), y=np.ravel(snr_iu1_group), color='blue')
    # plt.title('?? sessions \n TBC')
    plt.ylabel('SNR')
    plt.xlabel('FC')
    plt.savefig('/dhcp/fmri_anna_graham/GKgit/head_position/FINGER/finger_figures/fingerfigures_misc/SNR_covary_FC_kde.jpg')

