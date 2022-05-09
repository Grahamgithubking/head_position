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

    # print(allpid)
    # print(allsessid)
    # print()
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
    np.save('/dhcp/fmri_anna_graham/GKgit/fingerprinting/allpid.npy', allpid)
    print('The allpid list was saved as an array')
    np.save('/dhcp/fmri_anna_graham/GKgit/fingerprinting/allsessid.npy', allsessid)
    print('The allsessid dict was saved as an array')
    print()

def load_arrays():
    allpid = np.load('allpid.npy')
    allsessid = np.load('allsessid.npy', allow_pickle=True) #I had to allow_pickle=True so as to be able to load allsessid
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
    
    ##Two SWITCHES for choosing participants with either one (true/false), two sessions (false/true) or both(true/true):
    one_sessid = False
    two_sessid = True

    make_PID_SID_arrays(one_sessid, two_sessid)

    load_arrays()

    



