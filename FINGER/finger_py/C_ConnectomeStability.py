#### All original code is copyright of Graham King, Trinity College Institute of Neuroscience.

## Code to create figures/plots of connectome stability (within participant spearman correlation across the two sessions)
## Last edited Oct 2022

import pandas as pd
import numpy as np
import glob
import os
import nilearn.signal
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
from scipy import stats
import pingouin as pg



def get_two_session_subjects():
    # Identify all dHCP subjects with two sessions

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

    return allpid, allsessid

def calc_fwd(dhcp_root, allsessid):
    ## Calculating mean_fwd for each session

    allsess_mean=[]
    allsubject_mean=[]
    allsubject_max=[]

    for pidind, (pid, sessions) in enumerate(allsessid.items()):
        temp_mean=[]
        # for pid, sessions in allsessid.items():
        print("Working on participant %s"%pid)
        for sidind, sid in enumerate(sessions):
        # for sid in sessions:
            print("Working on session %s"%sid)

            # Reading motion.tsv file           
            fwd_path = os.path.join(dhcp_root, 'sub-' + str(pid), 'ses-' + str(sid), 'func')
            df_fwd = pd.read_csv(os.path.join(fwd_path, 'sub-' + str(pid) + '_ses-' + str(sid) + "_motion.tsv"), sep='\t')

            ## Calculating the mean of framewise_displacement:
            sess_mean=df_fwd['framewise_displacement'].mean(axis=0, skipna=True)
            print('The mean fwd for this session is:', sess_mean)
            temp_mean.append(sess_mean)
            allsess_mean.append(sess_mean)
        ##Calculating the mean fwd between the 2 sessions for each of 44 participants:
        subject_mean=np.mean(temp_mean)
        allsubject_mean.append(subject_mean)
        ##Calculating the max mean_fwd of the 2 sessions:
        subject_max=np.max(temp_mean)
        allsubject_max.append(subject_max)   

    print('the length of the allsess_mean list is:', len(allsess_mean))
    print('the length of the allsubject_mean list is:', len(allsubject_mean))
    print('the lenght of the allsubject_max list is:', len(allsubject_max))

    return allsubject_mean, allsubject_max

def create_dataframe(allpid, allsubject_mean, allsubject_max):
    # Creating dataframe

    withinfc=np.load('/dhcp/fmri_anna_graham/GKgit/finger_npy/44withinfc.npy') # These are the 44 Connectome Stability values (within participant spearman correlation across sessions)
    print('the length of the withinfc array is:', len(withinfc))
    print()

    snrmismean=np.load('/dhcp/fmri_anna_graham/GKgit/snr_npy/88_snr_416_erode1_mean.npy') # These are the 44 mean predicted SNR values (mean across 2 sessions and 387 nodes)
    print('the length of the snrmismean array is:', len(snrmismean))
    print()

    acrossmean= np.load('/dhcp/fmri_anna_graham/GKgit/finger_npy/44_mFC.npy')
    print('the length of the 44btwfc array is:', len(acrossmean))
    print()

    delta_mean = np.load('/dhcp/fmri_anna_graham/GKgit/finger_npy/44_deltaFC.npy')

    datacolumns = {'Participants':(allpid), 'Mean FWD':(allsubject_mean), 'Max FWD':(allsubject_max), 'Within Subj fc':(withinfc), 'Snrmismean':(snrmismean), '48_mFC':(acrossmean), '48_deltaFC':(delta_mean)}
    df = pd.DataFrame(datacolumns)
    # df.set_index('Participants')
    pd.set_option('display.max_rows', df.shape[0]+1)
    print(df)
    print()

    return withinfc, df


def adding_age_data(df):
    ##Adding in participant age data:

    df_age=pd.read_csv('/dhcp/fmri_anna_graham/GKgit/head_position/FINGER/finger_data/48subjbirth_noRS.csv')
    print('The imported birth dataframe is:')
    print(df_age)
    print()

    df_merge= pd.merge(left=df, right=df_age, left_on='Participants', right_on='id')
    pd.set_option('display.max_rows', df_merge.shape[0]+1)
    print('The combined dataframes are:')
    print(df_merge)
    print()

    #Saving dataframe to .csv file:
    df_merge.to_csv('/dhcp/fmri_anna_graham/GKgit/head_position/FINGER/finger_data/df_merge_withinfc.csv')
    print('The merged dataframe has been saved to .csv file')
    print()


    agelist=df_merge['birthgest']
    scan1list=df_merge['scan1age']
    scan2list=df_merge['scan2age']
    scanint=df_merge['scanint']

    return df_merge, agelist, scan1list, scan2list, scanint


def gettingplots(allsubject_mean, allsubject_max, withinfc, agelist, scan1list, scan2list, scanint, df_merge):
    # Scatter plots of Connectome Stability values versus others:

    plt.figure(1)
    norm=Normalize()
    norm.autoscale(scan1list)
    plt.scatter(x=(agelist), y=(withinfc), c=scan1list, norm = norm)
    cbar1=plt.colorbar()
    cbar1.ax.set_ylabel('Session 1 Age (weeks)')
    plt.title('Connectome Stability vs Birth Gestation')
    plt.xlabel('Birth Gestation (PMA weeks)')
    plt.ylabel('Connectome Stability')
    plt.savefig('/dhcp/fmri_anna_graham/GKgit/head_position/FINGER/finger_figures/fingerfigures_misc/birthgest_withinfp_sess1age.jpg')

    plt.figure(2)
    norm=Normalize()
    norm.autoscale(allsubject_mean)
    plt.scatter(x=(agelist), y=(withinfc), c=allsubject_mean, norm = norm)
    cbar2=plt.colorbar()
    cbar2.ax.set_ylabel('Mean framewise displacement')
    plt.title('Connectome Stability vs Birth Gestation')
    plt.xlabel('Birth Gestation (PMA weeks)')
    plt.ylabel('Connectome Stability')
    plt.savefig('/dhcp/fmri_anna_graham/GKgit/head_position/FINGER/finger_figures/fingerfigures_misc/birthgest_withinfc_meanfwd.jpg')

    plt.figure(3)
    plt.scatter(x=(scanint), y=(withinfc), c='grey')
    plt.xticks(range(0,16,2))
    plt.xlim([0, 16])
    plt.title('a)', fontsize=15, fontweight="bold")
    plt.xlabel('Inter-session interval (weeks)', fontsize=13)
    plt.ylabel('Connectome Stability', fontsize=13)
    plt.legend(['r=-0.47, *p<0.001'], loc='upper right', fontsize=13)
    plt.savefig('/dhcp/fmri_anna_graham/GKgit/head_position/FINGER/finger_figures/acrosssess_figures/CS_vs_scanint.jpg')

    plt.figure(4)
    plt.scatter(x=(allsubject_mean), y=(withinfc), c='grey')
    # plt.xticks(range(0,16,2))
    # plt.xlim([0, 16])
    # plt.title('Connectome Stability vs Mean FWD', fontsize=15, fontweight="bold")
    plt.xlabel('Mean FWD', fontsize=13)
    plt.ylabel('Connectome Stability', fontsize=13)
    plt.legend(['r=-0.5, *p<0.0005'], loc='upper right', fontsize=13)
    plt.savefig('/dhcp/fmri_anna_graham/GKgit/head_position/FINGER/finger_figures/acrosssess_figures/CS_vs_mFWD.jpg')


#######   Statistics, Partial Correlation coefficients_pearson with a single control:

    pc_one = pg.corr(scanint, withinfc, alternative='two-sided', method='pearson')
    print('The pearson correlation coefficient for Scan Interval and Within-Subject Connectome is:')
    print(pc_one)
    print()

    ppc_1a = pg.partial_corr(data=df_merge, x=('scanint'), y=('Within Subj fc'), covar=['Snrmismean'], alternative='two-sided', method='pearson')
    print('The partial correlation between Scan Interval and Within Subj fc, controlling for mean_SNRgroup, is:')
    print(ppc_1a)
    print()

    ppc_1b = pg.partial_corr(data=df_merge, x=('scanint'), y=('Within Subj fc'), covar=['48_deltaFC'], alternative='two-sided', method='pearson')
    print('The partial correlation between Scan Interval and Within Subj fc, controlling for 48_deltaFC, is:')
    print(ppc_1b)
    print()

    pc_two = pg.corr(allsubject_mean, withinfc, alternative='two-sided', method='pearson')
    print('The pearson correlation coefficient for Mean FWD and Within-Subject Connectome is:')
    print(pc_two)
    
    ppc_two = pg.partial_corr(data=df_merge, x=('Mean FWD'), y=('Within Subj fc'), covar=['Snrmismean'], alternative='two-sided', method='pearson')
    print('The partial correlation between Mean FWD and Within Subj fc, controlling for mean_SNRgroup, is:')
    print(ppc_two)
    print()

    pc_three = pg.corr(agelist, withinfc, alternative='two-sided', method='pearson')
    print('The pearson correlation coefficient for Birth Gestation and Within-Subject Connectome is:')
    print(pc_three)
    print()


    # pc_four = pg.partial_corr(data=df_merge, x=('birthgest'), y=('Within Subj fc'), covar=('Max FWD'), tail='two-sided', method='pearson')
    # print('The partial correlation between Birth Gestation and Within Subj fc, controlling for Max FWD, is:')
    # print(pc_four)
    # print()

    # pc_five = pg.partial_corr(data=df_merge, x=('Mean FWD'), y=('Within Subj fc'), covar=('scan1age'), tail='two-sided', method='pearson')
    # print('The partial correlation between Mean FWD and Within Subj fc, controlling for Session 1 Gestation, is:')
    # print(pc_five)
    # print()
    
    # pc_six = pg.partial_corr(data=df_merge, x=('Mean FWD'), y=('Within Subj fc'), covar=('scan2age'), tail='two-sided', method='pearson')
    # print('The partial correlation between Mean FWD and Within Subj fc, controlling for Session 2 Gestation, is:')
    # print(pc_six)
    # print()

    # pc_seven = pg.partial_corr(data=df_merge, x=('Mean FWD'), y=('Within Subj fc'), covar=('birthgest'), tail='two-sided', method='pearson')
    # print('The partial correlation between Mean FWD and Within Subj fc, controlling for Birth Gestation, is:')
    # print(pc_seven)
    # print()

    # pc_eight = pg.partial_corr(data=df_merge, x=('scanint'), y=('Within Subj fc'), covar=['scan1age'], tail='two-sided', method='pearson')
    # print('The partial correlation between Scan Interval and Within Subj fc, controlling for scan1age, is:')
    # print(pc_eight)
    # print()

    # pc_nine = pg.partial_corr(data=df_merge, x=('Mean FWD'), y=('Within Subj fc'), covar=['scan1age'], tail='two-sided', method='pearson')
    # print('The partial correlation between Mean FWD and Within Subj fc, controlling for scan1age, is:')
    # print(pc_nine)
    # print()

    # pc_ten = pg.partial_corr(data=df_merge, x=('Snrmismean'), y=('Within Subj fc'), covar=['scan1age'], tail='two-sided', method='pearson')
    # print('The partial correlation between Mean Predicted SNR and Within Subj fc, controlling for scan1age, is:')
    # print(pc_ten)
    # print()

if __name__ == '__main__':
    
    dhcp_root='/dhcp/dhcp_fmri_pipeline'
    
    allpid, allsessid= get_two_session_subjects()
   
    allsubject_mean, allsubject_max= calc_fwd(dhcp_root, allsessid)

    withinfc, df = create_dataframe(allpid, allsubject_mean, allsubject_max)

    df_merge, agelist, scan1list, scan2list, scanint = adding_age_data(df)
    
    gettingplots(allsubject_mean, allsubject_max, withinfc, agelist, scan1list, scan2list, scanint, df_merge)

    