#### All original code is copyright of Graham King, Trinity College Institute of Neuroscience.

### Code to creat various numpy arrays of specific participant ID/session ID from participants with either one, two sessions, or both:
###NOTE: The Switches to select participants with either one or two sessions
###NOTE: The first N number of participants in the first dataframe

import pandas as pd
import numpy as np
import glob
import os

def make_PID_SID_arrays(one_sessid, two_sessid):

    dhcp_root='/dhcp/dhcp_fmri_pipeline'

    # Load list of subjects
    df_all=pd.read_csv(os.path.join(dhcp_root,'participants.tsv'), delimiter='\t')
    # df_all=df_all.head(20) #Testing first N number of participants

    # Load list of sessions
    allpid=[]
    print(type(allpid))
    allsessid={}
    print(type(allsessid))
    
    for pid in df_all['participant_id']: 
        df_sess=pd.read_csv(os.path.join(dhcp_root, 'sub-' + pid, 'sub-' + pid + "_sessions.tsv" ), delimiter='\t')
        if one_sessid:
            if (len(df_sess))<2:
            #Choosing partcipants with only one sessid
                allpid.append(pid)
                allsessid[pid]=[]
                for sid in df_sess['session_id']:
                    allsessid[pid].append(sid)
        if two_sessid:
            if (len(df_sess))>1:
            #Choosing partcipants with two sessid
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

    print('The length of allpid is:')
    print(len(allpid))
    maxnsess=np.max([len(x)for x in allsessid.values()])
    print('The length of maxnsess per participant is:')
    print(maxnsess)
    print()

    allpid = np.array(allpid)
    allsessid = np.array(list(allsessid.items()), dtype=dict) #allsessid is now an array
    print('The allsessid dict has been converted to array of type:')
    print(type(allsessid))
    print()

    if one_sessid:
        np.save('/dhcp/fmri_anna_graham/GKgit/finger_npy/416allpid.npy', allpid)
        np.save('/dhcp/fmri_anna_graham/GKgit/finger_npy/416allsessid.npy', allsessid)
    if two_sessid:
        np.save('/dhcp/fmri_anna_graham/GKgit/finger_npy/44allpid.npy', allpid)
        np.save('/dhcp/fmri_anna_graham/GKgit/finger_npy/88allsessid.npy', allsessid)


    print('The allpid list was saved as an array')
    print('The allsessid dict was saved as an array')
    print()

def load_arrays(one_sessid, two_sessid):
    
    if one_sessid:
        allpid = np.load('/dhcp/fmri_anna_graham/GKgit/finger_npy/416allpid.npy')
        allsessid = np.load('/dhcp/fmri_anna_graham/GKgit/finger_npy/416allsessid.npy', allow_pickle=True)
    if two_sessid:
        allpid = np.load('/dhcp/fmri_anna_graham/GKgit/finger_npy/44allpid.npy')
        allsessid = np.load('/dhcp/fmri_anna_graham/GKgit/finger_npy/88allsessid.npy', allow_pickle=True)
    print('The two npy array files were loaded')
    print()

    print('This is allpid:')
    print(allpid)
    print(type(allpid))
    print('This is allsessid:')
    print(allsessid) 
    print(type(allsessid)) 

    allsessid = dict(allsessid)
    # allsessid = dict(np.ndenumerate(allsessid))
    print('Reconverted allsessid to type dict:')
    # print(allsessid)
    print(type(allsessid))
                
    
if __name__ == '__main__':
    
    ##Two SWITCHES for choosing participants with either one (true/false), two sessions (false/true):
    one_sessid = False
    two_sessid = True

    make_PID_SID_arrays(one_sessid, two_sessid)

    load_arrays(one_sessid, two_sessid)

    
