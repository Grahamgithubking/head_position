#### All original code is copyright of Graham King, Trinity College Institute of Neuroscience.
### Code to calculate the Spearman second order correlation
### Code to calculate fingerprint identification success
### Code to calculate mean within vs between subject correlation
### Fishers exact permutation test to determine significance
### Last edited by GK on June 14 2021

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
    print(allpid)
    print()
    return allpid, allsessid

def session2session(allpid):
    
    # loading allfc
    allfc = np.load('/dhcp/fmri_anna_graham/GKgit/finger_npy/fc_96.npy')
    sz=allfc.shape
    print('starting allfc shape is:')
    print(sz)
    nsubj=sz[0]
    nsess=sz[1]
    nroi=sz[2]

    # Reshaping in full 387x387 matrix:
    allfc_reshaped=np.reshape(allfc, (nsubj * nsess, nroi * nroi))
    print('The shape of allfc_reshaped is:')
    print(allfc_reshaped.shape)
    
    # # Reshaping in Upper Triangle:
    # iu1=np.triu_indices(nroi, k=1)
    # allfc_reshaped=np.reshape(allfc, (nsubj*nsess, nroi, nroi))
    # print('The shape of allfc_reshaped is:')
    # print(allfc_reshaped.shape)
    # allfc_iu1=np.array([x[iu1]for x in allfc_reshaped]) #for each session, selecting the upper triangle of the 387x387 matrix
    # print('The shape of allfc_iu1 is:')
    # print(allfc_iu1.shape)
    # allfc_reshaped = allfc_iu1

    ## Spearman second order correlation across 96 sessions:
    sessim, pval=stats.spearmanr(allfc_reshaped.T)
    print('the shape of the sessim matrix is:')  # this is a 96x96 matrix
    print(sessim.shape)
    print()


    ########### Fingerprint Identification analysis ######################################:

    #Find all the points comparing session 0 to session 1
    comparesess=sessim[0::2,1::2]
    print('The shape of comparesess is:')
    print(comparesess.shape)
    print('Sliced at [42:,42:] gives:')
    print()
    print(comparesess[42:,42:])
    print()


    ## Saving the within participant array:
    within=np.diag(comparesess)
    print(within)
    print('The shape of the within array is:')
    print(within.shape)
    print()
    np.save('/dhcp/fmri_anna_graham/GKgit/finger_npy/48withinfc.npy', within)

    ## Saving the across particpant array:
    btwfc_all=[]
    for i in range (nsubj):
        i_row = comparesess[i]
        btw = np.delete(i_row, i)
        btwfc = np.mean(btw)
        btwfc_all.append(btwfc)
    btwfc_array=np.array(btwfc_all)
    print(btwfc_array)
    print('The shape of the btwfc_all array is:')
    print(btwfc_array.shape)
    print()
    np.save('/dhcp/fmri_anna_graham/GKgit/finger_npy/btwfc.npy', btwfc_array)
    
    ## Saving the unique fingerprint fp_i (within_fc minus between_fc, for each i th subject):
    fp_all=[]
    for i in range (nsubj):
        fc_i = comparesess[i,i]
        i_row = comparesess[i]
        less_fc_i = np.delete(i_row, i)
        fp_i = fc_i - np.mean(less_fc_i)
        fp_all.append(fp_i)
    fp_array=np.array(fp_all)
    print(fp_all)
    print('the shape of the fp_all array is:')
    print(fp_array.shape)
    print()
    np.save('/dhcp/fmri_anna_graham/GKgit/finger_npy/48fp.npy', fp_array)
    
    

    plt.figure(figsize=(12,8))
    c = plt.imshow(comparesess, cmap='Greens') #Edit color
    plt.colorbar(c)
    plt.title('a)', fontsize=20, fontweight="bold")
    plt.xlabel('Participants 1 - 48', fontsize=18)
    plt.ylabel('Participants 1 - 48', fontsize=18)
    plt.savefig('/dhcp/fmri_anna_graham/GKgit/head_position/FINGER/finger_figures/acrosssess_figures/RSM_acrosssess.jpg')

    #Sort columns according to how well each subject of the session0 matches each of the session1
    comparesess_sorted=np.argsort(comparesess, axis=0)
    print('The shape of comparesess_sorted is:')
    print(comparesess_sorted.shape)
    print('Sliced at [42:,42:] gives:')
    print()
    print(comparesess_sorted[42:,42:])
    print()

    # Identifying whether withinFC is the highest correlation:
    match = np.diag(comparesess) == np.max(comparesess, axis=0)
    for ind, subj in enumerate(allpid):
        print('%s\t%s'%(subj, match[ind]))
    print()

    #Find out the rank of the true match (in column i, where subject i ended up in the sorted ranking)
    ### Error with code here????????
    rankofmatch=np.where(((comparesess_sorted - np.arange(nsubj))==0))[0]
    for ind, subj in enumerate(allpid):
        print('%s\t%d'%(subj, rankofmatch[ind]))
    print()


    bins=np.arange(0,49,1)-0.5
    plt.style.use('ggplot')
    plt.figure(figsize=(16,12))
    plt.hist(rankofmatch, bins, rwidth=0.5, color='g')
    plt.xticks(range(0,48,1))
    plt.xlim([-1, 49])
    plt.yticks(range(0,10,1))
    plt.title('Rank of the within-subject correlation across sessions', fontsize=20, pad=5, loc='center')
    plt.xlabel('Rank (lowest = 0 to highest = 47)', fontsize=15)
    plt.ylabel('No. Subjects', fontsize=15)
    plt.savefig('/dhcp/fmri_anna_graham/GKgit/fingerprinting/FINGER/finger_figures/acrosssess_figures/ID_sessions.png')
    plt.show()

    ## Calculating the difference in correlation meanwithin vs meanbetween:
    iu1 = np.triu_indices(nsubj, k=1)
    within = np.diag(comparesess)
    between = comparesess[iu1]

    N_within=len(within)
    mn_within=np.mean(within)
    se_within=np.std(within)/np.sqrt(N_within)
    print("Within subject comparison {} values {} +/- {}".format(N_within, mn_within, se_within))
    print()

    N_across=len(between)
    mn_across=np.mean(between)
    se_across=np.std(between)/np.sqrt(N_across)
    print("Across subject comparison {} values {} +/- {}".format(N_across, mn_across, se_across))
    print()

    mean_within_minus_between = np.mean(within) - np.mean(between)
    print(f'Second Level Correlation - Mean within subject minus Mean between subject:{mean_within_minus_between}')
    print()

    return comparesess, mean_within_minus_between, nsubj, iu1
    
def permutation_test(comparesess, mean_within_minus_between, nsubj, iu1):
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

    plt.figure(figsize=(16,12))
    bins=np.arange(0,49,1)-0.5
    plt.hist(alltopshuffle, bins, rwidth=0.5, color='b')
    plt.xticks(range(0,48,1))
    plt.xlim([-1, 49])
    plt.yticks(range(0,10000,50))
    plt.ylim([0, 500])
    plt.title('Fisher 10k - Rank', fontsize=20)
    plt.xlabel('No. of Top ranks', fontsize=15)
    plt.ylabel('No. Permutations', fontsize=15)
    plt.savefig('/dhcp/fmri_anna_graham/GKgit/fingerprinting/FINGER/finger_figures/acrosssess_figures/Rank_Fisher.png') #edit name here
    plt.show()


    plt.figure(figsize=(16,12))
    bins=np.arange(0,0.1,0.005)
    plt.hist(shuffle_mean_within_minus_between, bins, rwidth=0.5, color='b')
    plt.yticks(range(0,10000,50))
    plt.ylim([0, 3000])
    plt.title('Fisher 10k - DeltaR', fontsize=20)
    plt.xlabel('Delta r values', fontsize=15)
    plt.ylabel('No. Permutations', fontsize=15)
    plt.savefig('/dhcp/fmri_anna_graham/GKgit/fingerprinting/FINGER/finger_figures/acrosssess_figures/DeltaR_Fisher.png') #edit name here
    plt.show()

    within_between_perm = np.mean(shuffle_mean_within_minus_between>=mean_within_minus_between) ##No. of trues divided by 10k
    print(f'Shuffled - Mean within subject minus Mean between subject, p<{within_between_perm}')
    

if __name__ == '__main__':

    timecourse_pth = '/dhcp/fmri_anna_graham/timecourses/'
    allpid, allsessid = get_two_session_subjects()
    comparesess, mean_within_minus_between, nsubj, iu1 = session2session(allpid)
    permutation_test(comparesess, mean_within_minus_between,nsubj, iu1)




