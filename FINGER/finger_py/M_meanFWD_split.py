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

    all_fwd_split=np.empty((44,2,2))

    for pidind, (pid, sessions) in enumerate(allsessid.items()):
        # for pid, sessions in allsessid.items():
        print("Working on participant %s"%pid)
        for sidind, sid in enumerate(sessions):
            # for sid in sessions:
            print("Working on session %s"%sid)

            # Reading motion.tsv file           
            fwd_path = os.path.join(dhcp_root, 'sub-' + str(pid), 'ses-' + str(sid), 'func')
            df_fwd = pd.read_csv(os.path.join(fwd_path, 'sub-' + str(pid) + '_ses-' + str(sid) + "_motion.tsv"), sep='\t')

            fwd=df_fwd['framewise_displacement']
            # fwd=list(fwd)
            fwd=np.array(fwd)
            print(f"The length of all fwd for session is: {len(fwd)}")
            print(f"The shape of all fwd for session is: {fwd.shape}")
            print()

            # Remove 60sec 154volumes
            num_range = range(1073, 1227) #zero based
            print(len(num_range))
            crop = list(num_range) 
            fwd_clipped = np.delete(fwd, crop)
            print('The length and shape of fwd_clipped is:')
            print(len(fwd_clipped))
            print(fwd_clipped.shape)

            #split timecourse into two
            fwd_split = np.split(fwd_clipped, 2)
            split_A = fwd_split[0]
            split_B = fwd_split[1]

            print('This is the shape of split_A:')
            print(split_A.shape)
            print()

            split_A_mean=np.mean(split_A)
            split_B_mean=np.mean(split_B)
            print(f"The mean of split_A is: {split_A_mean}")
            print(f"The mean of split_B is: {split_B_mean}")
            print()

            all_fwd_split[pidind,sidind,0]=split_A_mean
            all_fwd_split[pidind,sidind,1]=split_B_mean

    np.save('/dhcp/fmri_anna_graham/GKgit/finger_npy/176_mFWD.npy', all_fwd_split)


if __name__ == '__main__':
    
    dhcp_root='/dhcp/dhcp_fmri_pipeline'

    allpid, allsessid= get_two_session_subjects()

    calc_fwd(dhcp_root, allsessid)
