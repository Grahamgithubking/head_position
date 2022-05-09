#### All original code is copyright of Graham King, Trinity College Institute of Neuroscience.
## Code to create Split Segments of each of 96 session for 48 Participants

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
 

def split_timecourses(timecourse_pth, split_pth, allsessid):
    # Aiming to split each session into half!

    for pidind, (pid, sessions) in enumerate(allsessid.items()):
    # for pid, sessions in allsessid.items():
        print("Working on participant %s"%pid)
        for sidind, sid in enumerate(sessions):
        # for sid in sessions:
            print("Working on session %s"%sid)

            # Load up timecourse
            filename = os.path.join(timecourse_pth, 'sub-' + str(pid), 'ses-' + str(sid), 'sub-' + str(pid) + '_ses-' + str(sid) + "_task-rest_bold.txt")
            print("Loading %s"%filename)
            timecourse = np.loadtxt(filename)
            print(timecourse.shape)

            # Remove 60sec 154volumes
            num_range = range(1073, 1227) #zero based
            print(len(num_range))
            crop = list(num_range) 
            timecourse_clipped = np.delete(timecourse, crop, axis=0) #axis0 is time
            print('This is timecourse_clipped:')
            print(timecourse_clipped.shape)
            print()

            #split timecourse into two
            timecourse_split = np.vsplit(timecourse_clipped, 2)
            split_A = timecourse_split[0]
            split_B = timecourse_split[1]

            print('This is the size of split_A:')
            print(split_A.shape)
            print()

            np.save(split_pth + str(pid) + str(sid) + '_A', split_A)
            np.save(split_pth + str(pid) + str(sid) + '_B', split_B)


if __name__ == '__main__':
    allpid, allsessid= get_two_session_subjects()

    timecourse_pth = '/dhcp/fmri_anna_graham/timecourses/'
    split_pth = '/dhcp/fmri_anna_graham/timecourses_split/'

    split_timecourses(timecourse_pth, split_pth, allsessid)
