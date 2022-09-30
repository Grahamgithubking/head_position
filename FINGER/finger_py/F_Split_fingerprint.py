#### All original code is copyright of Graham King, Trinity College Institute of Neuroscience.
### Code to calculate the Spearman second order correlation
### Code to calculate fingerprint identification success
### Code to calculate mean within vs between subject correlation
### Fishers exact permutation test to determine significance
### Last edited by GK on April 7th 2021

import pandas as pd
import numpy as np
import glob
import os
import nilearn.signal
import matplotlib.pyplot as plt
from scipy import stats



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

    return allpid, allsessid


def session2session(allpid, oldersession):
    # loading split segment first order pearson correlations as allfc
    allfc = np.load('/dhcp/fmri_anna_graham/GKgit/finger_npy/192fc.npy')
    sz=allfc.shape
    print('starting allfc shape is:')
    print(sz)
    nsubj=sz[0]
    nsess=sz[1]
    nsplit=sz[2]
    nroi=sz[3]

    ## Reshaping 
    allfc_reshaped=np.reshape(allfc, (nsubj * nsess * nsplit, nroi * nroi))
    print('The shape of allfc_reshaped is:')
    print(allfc_reshaped.shape)


    ## Spearman second order correlation across 192 split segments:
    sessim, pval=stats.spearmanr(allfc_reshaped.T)
    print('the shape of the sessim matrix is:') # this is a 192x192 matrix
    print(sessim.shape)
    print()



    ########### Fingerprint Identification analysis ######################################:

    #Find all the points comparing split segment A to split segment B
    if oldersession:
        comparesess=sessim[2::4,3::4] #this chooses the older session
        print('the shape of older comparesess is:')
        print(comparesess.shape)
    else:
        comparesess=sessim[0::4,1::4] #this chooses the younger session
        print('the shape of younger comparesess is:')
        print(comparesess.shape)

    if oldersession:
        plt.figure(figsize=(12,8))
        c = plt.imshow(comparesess, cmap='Reds') # Edit color here
        plt.colorbar(c)
        plt.title('c)', fontsize=20, fontweight="bold") # Edit title
        plt.xlabel('Participants 1 - 48', fontsize=18)
        plt.ylabel('Participants 1 - 48', fontsize=18)
        plt.savefig('/dhcp/fmri_anna_graham/GKgit/head_position/FINGER/finger_figures/fc_figures/RSM_split_older.jpg') #edit older/younger here!
    else:
        plt.figure(figsize=(12,8))
        c = plt.imshow(comparesess, cmap='Blues') # Edit color here
        plt.colorbar(c)
        plt.title('b)', fontsize=20, fontweight="bold") # Edit title
        plt.xlabel('Participants 1 - 48', fontsize=18)
        plt.ylabel('Participants 1 - 48', fontsize=18)
        plt.savefig('/dhcp/fmri_anna_graham/GKgit/head_position/FINGER/finger_figures/fc_figures/RSM_split_younger.jpg') #edit older/younger here!

    # #Sort columns according to how well each subject of split A matches each of split B
    # comparesess_sorted=np.argsort(comparesess, axis=0)
    # #Find out the rank of the true match (in column i, where subject i ended up in the sorted ranking)
    # rankofmatch=np.where((comparesess_sorted - np.arange(nsubj))==0)[0]
    # for ind, subj in enumerate(allpid):
    #     print('%s\t%d'%(subj, rankofmatch[ind]))
    # print()

    # Identifying whether withinFC is the highest correlation:
    match = np.diag(comparesess) == np.max(comparesess, axis=0)
    for ind, subj in enumerate(allpid):
        print('%s\t%s'%(subj, match[ind]))
    print()


    ## Calculating the difference in correlation meanwithin vs meanbetween:
    iu1 = np.triu_indices(nsubj, k=1)
    within = np.diag(comparesess)
    between = comparesess[iu1]
    if oldersession:
        np.save('/dhcp/fmri_anna_graham/GKgit/finger_npy/48withinfc_splitsterm.npy', within)
    else:
        np.save('/dhcp/fmri_anna_graham/GKgit/finger_npy/48withinfc_splitspreterm.npy', within)

    mean_within = np.mean(within)
    print(f'Second level correlation - Mean within subject:{mean_within}')
    mean_between = np.mean(between)
    print(f'Second level correlation - Mean between subject:{mean_between}')
    mean_within_minus_between = mean_within - mean_between
    print(f'Second Level Correlation - Mean within subject minus between subject:{mean_within_minus_between}')
    print()

    return comparesess, mean_within_minus_between, nsubj, iu1

def permutation_test(comparesess, mean_within_minus_between, nsubj, iu1, oldersession):
    # Fishers Exact Test with 10,000 permutations:    
    alltopshuffle=[]
    shuffle_mean_within_minus_between=[]
    shufflerows=comparesess
    
    for perm in range(9999):
        np.random.shuffle(shufflerows)
        shufflerows_sorted=np.argsort(shufflerows, axis=0)
        rankofshuffle=np.where((shufflerows_sorted - np.arange(nsubj))==0)[0]

        ##permutations for toprank:
        topshuffle = np.where(rankofshuffle>=47)[0]
        lengthtopshuffle = len(topshuffle)
        alltopshuffle.append(lengthtopshuffle)

        ##permutations for within_versus_between:
        within = np.diag(shufflerows) #edited from comparesess
        between = shufflerows[iu1]    #edited from comparesess
        shuffle_mean_within_minus_between.append(np.mean(within) - np.mean(between))

    if oldersession:
        plt.figure(figsize=(16,12))
        bins=np.arange(0,49,1)-0.5
        plt.hist(alltopshuffle, bins, rwidth=0.5, color='b')
        plt.xticks(range(0,48,1))
        plt.xlim([-1, 49])
        plt.yticks(range(0,10000,50))
        plt.ylim([0, 500])
        plt.title('Fisher 10k - Rank Older', fontsize=20)
        plt.xlabel('No. of Top ranks', fontsize=15)
        plt.ylabel('No. Permutations', fontsize=15)
        plt.savefig('/dhcp/fmri_anna_graham/GKgit/head_position/FINGER/finger_figures/fc_figures/Rank_Fisher_Older.png') #edit name here
        plt.show()
    else:
        plt.figure(figsize=(16,12))
        bins=np.arange(0,49,1)-0.5
        plt.hist(alltopshuffle, bins, rwidth=0.5, color='b')
        plt.xticks(range(0,48,1))
        plt.xlim([-1, 49])
        plt.yticks(range(0,10000,50))
        plt.ylim([0, 500])
        plt.title('Fisher 10k - Rank Younger', fontsize=20)
        plt.xlabel('No. of Top ranks', fontsize=15)
        plt.ylabel('No. Permutations', fontsize=15)
        plt.savefig('/dhcp/fmri_anna_graham/GKgit/head_position/FINGER/finger_figures/fc_figures/Rank_Fisher_Younger.png') #edit name here
        plt.show()

    if oldersession:
        plt.figure(figsize=(16,12))
        bins=np.arange(0,0.1,0.005)
        plt.hist(shuffle_mean_within_minus_between, bins, rwidth=0.5, color='b')
        plt.yticks(range(0,10000,50))
        plt.ylim([0, 3000])
        plt.title('Fisher 10k - Delta r Older', fontsize=20)
        plt.xlabel('Delta r values', fontsize=15)
        plt.ylabel('No. Permutations', fontsize=15)
        plt.savefig('/dhcp/fmri_anna_graham/GKgit/head_position/FINGER/finger_figures/fc_figures/DeltaR_Fisher_Older.png') #edit name here
        plt.show()
    else:
        plt.figure(figsize=(16,12))
        bins=np.arange(0,0.1,0.005)
        plt.hist(shuffle_mean_within_minus_between, bins, rwidth=0.5, color='b')
        plt.yticks(range(0,10000,50))
        plt.ylim([0, 3000])
        plt.title('Fisher 10k - Delta r Younger', fontsize=20)
        plt.xlabel('Delta r values', fontsize=15)
        plt.ylabel('No. Permutations', fontsize=15)
        plt.savefig('/dhcp/fmri_anna_graham/GKgit/head_position/FINGER/finger_figures/fc_figures/DeltaR_Fisher_Younger.png') #edit name here
        plt.show()

    within_between_perm = np.mean(shuffle_mean_within_minus_between>=mean_within_minus_between) ##No. of trues divided by 10k
    print(f'Shuffled - Mean within subject minus Mean between subject, p<{within_between_perm}')


if __name__ == '__main__':
    
    ##SWITCH for choosing younger vs older session:
    oldersession = True

    allpid, allsessid= get_two_session_subjects()
    comparesess, mean_within_minus_between, nsubj, iu1 = session2session(allpid, oldersession)
    permutation_test(comparesess, mean_within_minus_between, nsubj, iu1, oldersession)



