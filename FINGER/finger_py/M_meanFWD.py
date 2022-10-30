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
            print()
            temp_mean.append(sess_mean)
            allsess_mean.append(sess_mean)
        ##Calculating the mean fwd between the 2 sessions for each of 48 participants:
        subject_mean=np.mean(temp_mean)
        allsubject_mean.append(subject_mean)
        ##Calculating the max mean_fwd of the 2 sessions:
        subject_max=np.max(temp_mean)
        allsubject_max.append(subject_max)   

    print(allsubject_mean)
    print()

    np.save('/dhcp/fmri_anna_graham/GKgit/finger_npy/44_mFWD.npy', allsubject_mean)

if __name__ == '__main__':
    
    dhcp_root='/dhcp/dhcp_fmri_pipeline'
    
    allpid, allsessid= get_two_session_subjects()

    calc_fwd(dhcp_root, allsessid)
