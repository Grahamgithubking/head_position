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


def get_two_session_subjects():
    # Identify all dHCP subjects with two sessions

    dhcp_root='/dhcp/dhcp_fmri_pipeline'

    # Load list of subjects
    df_subj=pd.read_csv(os.path.join(dhcp_root,'participants.tsv'), delimiter='\t')

    # Load list of sessions
    allpid=[]
    allsessid={}

    for pid in df_subj['participant_id']:
        df_sess=pd.read_csv(os.path.join(dhcp_root, 'sub-' + pid, 'sub-' + pid + "_sessions.tsv" ), delimiter='\t')
        if (len(df_sess))>1:
            # More than one session
            allpid.append(pid)
            allsessid[pid]=[]
            for sid in df_sess['session_id']:
                allsessid[pid].append(sid)

    del_pid=['CC00191XX11', 'CC00518XX15', 'CC00672AN13', 'CC00770XX12']
    for pid in del_pid:
        allpid.remove(pid)
        del allsessid[pid]
    print(f"This is the updated allpid: \n {allpid}")
    print(f"This is the updated allsessid: \n {allsessid}")
    print()

    return allpid

def loadsnr(imageroot, one_full, two_full, two_split, snrcoil, snrtrue, oldersession, allpid):

    if one_full:
        if snrcoil:
            x='snrcoil'
            y='416'
            snr_all=np.load('/dhcp/fmri_anna_graham/GKgit/snr_npy/416_snr_416_erode1.npy') #This has shape [416,1,387]
        if snrtrue:
            x='snrtrue'
            y='416'
            snr_all=np.load('/dhcp/fmri_anna_graham/GKgit/snr_npy/416_snr_true.npy')

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

    if two_full:
        if snrcoil:
            x='snrcoil'
            y='88'
            snr_all=np.load('/dhcp/fmri_anna_graham/GKgit/snr_npy/88_snr_416_erode1.npy') #This has shape [44,2,387]
        if snrtrue:
            x='snrtrue'
            y='88'
            snr_all=np.load('/dhcp/fmri_anna_graham/GKgit/snr_npy/88_snr_true.npy') #This has shape [44,2,387]

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

    if two_split:
        if snrcoil:
            snr_all = np.load('/dhcp/fmri_anna_graham/GKgit/snr_npy/176hp_snr.npy')
            x = 'snrcoil'
            if oldersession:
                y='88split_older_'
            else:
                y='88split_younger_'
        if snrtrue:
            snr_all = np.load('/dhcp/fmri_anna_graham/GKgit/snr_npy/176_snrtrue.npy')
            x = 'snrtrue'
            if oldersession:
                y='88split_older_'
            else:
                y='88split_younger_'
                    

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



    #### I <did/did not> DO THIS:
    #### Multiplying the SNR across all edges for each ROI:
    if two_split:
        hp_all=nsubj*nsess*nhp
    else:
        hp_all=nsubj*nsess
    snrbysnr=np.zeros((hp_all,nroi,nroi))
    for ind in range(hp_all):
        snr=snr_all_reshaped[ind,:]
        snrbysnr[ind,:,:]=np.outer(snr, snr) #instead of using snr*snr.T
    print('The shape of snrbysnr is:')
    print(snrbysnr.shape)
    print()
    iu1=np.triu_indices(nroi, k=1) #selecting the upper triangle of a matrix
    allsnr_iu1=np.array([x[iu1]for x in snrbysnr]) #for each session, selecting the upper triangle of the 387x387 matrix
    print('The shape of allsnr_iu1 is:')
    print(allsnr_iu1.shape)
    snr_all_reshaped=allsnr_iu1


    snrim, pval=stats.spearmanr(snr_all_reshaped.T)
    print('The shape of the spearman correlation of snr_all is:')
    print(snrim.shape)
    print() 

    if two_full:
        comparesess=snrim[0::2,1::2]
        print('The shape of comparesess is:')
        print(comparesess.shape)
        print('Sliced at [40:,40:] gives:')
        print()
        print(comparesess[40:,40:])
        print()

        plt.figure(figsize=(12,8))
        c = plt.imshow(comparesess, cmap='Greens') #Edit color
        plt.colorbar(c)
        plt.title('a)', fontsize=20, fontweight="bold")
        plt.xlabel('Participants 1 - 44', fontsize=18)
        plt.ylabel('Participants 1 - 44', fontsize=18)
        

    if two_split:
        if oldersession:
            comparesess=snrim[2::4,3::4] #this chooses the older session
            print('the shape of older comparesess is:')
            print(comparesess.shape)
        else:
            comparesess=snrim[0::4,1::4] #this chooses the younger session
            print('the shape of younger comparesess is:')
            print(comparesess.shape)

        if oldersession:
            plt.figure(figsize=(12,8))
            c = plt.imshow(comparesess, cmap='Reds') # Edit color here
            plt.colorbar(c)
            plt.title('c)', fontsize=20, fontweight="bold") # Edit title
            plt.xlabel('Participants 1 - 44', fontsize=18)
            plt.ylabel('Participants 1 - 44', fontsize=18)
        else:
            plt.figure(figsize=(12,8))
            c = plt.imshow(comparesess, cmap='Blues') # Edit color here
            plt.colorbar(c)
            plt.title('b)', fontsize=20, fontweight="bold") # Edit title
            plt.xlabel('Participants 1 - 44', fontsize=18)
            plt.ylabel('Participants 1 - 44', fontsize=18)

    plt.savefig(str(imageroot) + str(y) + str(x) + '.jpg')

    # Identifying whether withinSNR is the highest correlation:
    match = np.diag(comparesess) == np.max(comparesess, axis=0)
    for ind, subj in enumerate(allpid):
        print('%s\t%s'%(subj, match[ind]))
    print()


    ## Calculating the difference in SNRcorrelation meanwithin vs meanbetween:
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

    imageroot='/dhcp/fmri_anna_graham/GKgit/head_position/FINGER/finger_figures/snr_figures/'

    ## SWITCHES:
    snrcoil=True
    snrtrue=False
    
    one_full=False
    two_full=True
    two_split=False
    oldersession=False

    allpid = get_two_session_subjects()
  
    loadsnr(imageroot, one_full, two_full, two_split, snrcoil, snrtrue, oldersession, allpid)

