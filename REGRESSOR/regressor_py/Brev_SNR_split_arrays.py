# ## All original code is copyright of Graham King, Trinity College Institute of Neuroscience. 
# #Date: 08/2022

# Inputting of <Predicted/True> SNR values into numpy arrays for two splits per session,
# for all sessions of the 48 particpants with both preterm and term sessions
# 
# allsnr is a 48x2x2x387 npy matrix

# TBC.............

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


def load_snr(hp_path, allsessid, snrtrue):
    # Loading SNRCoil or SNRTrue values for TWO splits in each of 96 sessions

    nrois=387
    nhps=2 # number of head positions/splits
    nsess=2
    nsubj=48
    allsnr=np.zeros((nsubj, nsess, nhps, nrois))
    print('the shape of allsnr zeros array is:')
    print(allsnr.shape)
    print()

    if snrtrue:
        hps = ['_SNRtruesplita', '_SNRtruesplitb']
    else:
        hps = ['_SNR_HP537', '_SNR_HP1763']

    for pidind, (pid, sessions) in enumerate(allsessid.items()):
    # for pid, sessions in allsessid.items():
        print("Working on participant %s"%pid)
        for sidind, sid in enumerate(sessions):
        # for sid in sessions:
            print("Working on session %s"%sid)
            for hpind, hp in enumerate(hps):
            #for hp in hps:
                print("Working on hp %s"%hp)

                #Load up SNR values for each hp:
                filename = os.path.join(hp_path, 'sub-' + str(pid), 'ses-' + str(sid), 'sub-' + str(pid) + '_ses-' + str(sid) + str(hp) + ".txt")
                print("Loading %s"%filename)
                print()
                snr = np.loadtxt(filename)

                # Remove bad columns
                rois2reject_onebased=[2, 3, 4, 6, 113, 116, 139, 203, 204, 205, 314, 319, 347]
                rois2reject_zerobased=[x-1 for x in rois2reject_onebased]
                snr_trimmed = np.delete(snr, rois2reject_zerobased)
                print('The shape of the snr_trimmed array is:')
                print(snr_trimmed.shape)
                print()  

                #Matrix for all 192 head position SNR values:
                allsnr[pidind, sidind, hpind, :]=snr_trimmed

    print('The shape of the 192_split_snr array is:')
    print(allsnr.shape)
    print()

    # Saving output:
    if snrtrue:
        np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/192_snrtrue.npy', allsnr)
    else:
        np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/192hp_snr.npy', allsnr)
    
    

if __name__ == '__main__':
    
    ## Switch to select snrcoil (false) versus snrtrue (true):
    snrtrue=False
    
    allpid, allsessid= get_two_session_subjects()

    hp_path = '/dhcp/scanner_positioning/rois/'
    
    load_snr(hp_path, allsessid, snrtrue)
    