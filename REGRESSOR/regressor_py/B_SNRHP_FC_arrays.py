#### All original code is copyright of Graham King, Trinity College Institute of Neuroscience.

###Step1: Loading predicted SNR-Headposition txt files and saving as numpy array
###Step2: Loading BOLD fmri txt files, removing QC ROIs, filtering, determining the functional connectivity (first order pearson), and saving as numpy array
###NOTE: Switches below to Turn on/off either step
###NOTE: Change name of saved output files!!!!

import pandas as pd
import numpy as np
import glob
import os
import nilearn.signal
import matplotlib.pyplot as plt
from scipy import stats
import pingouin as pg
import statsmodels.api as sm

def load_arrays():
    allpid = np.load('/dhcp/fmri_anna_graham/GKgit/finger_npy/allpid.npy')
    allsessid = np.load('/dhcp/fmri_anna_graham/GKgit/finger_npy/allsessid.npy', allow_pickle=True) #I had to allow_pickle=True so as to be able to load allsessid
    print('The two npy array files were loaded')
    print()

    print('This is allpid:')
    print(allpid)
    print(type(allpid))
 
    allsessid = dict(allsessid)
    print('Reconverted allsessid to type dict:')
    print(allsessid)
    print(type(allsessid))
    
    return allpid, allsessid


def loading_SNRHP(run_snr, snr_path, allpid, allsessid):

    maxnsess=np.max([len(x)for x in allsessid.values()])
    snr_all=np.zeros((len(allsessid), maxnsess, 387))
    print('The original shape of snr_all zeros is:')
    print(snr_all.shape)

    if run_snr:
        for pidind, (pid, sessions) in enumerate(allsessid.items()):
        # for pid, sessions in allsessid.items():
            print("Working on participant %s"%pid)
            for sidind, sid in enumerate(sessions):
            # for sid in sessions:
                print("Working on session %s"%sid)

                filename = os.path.join(snr_path, 'sub-' + str(pid), 'ses-' + str(sid), 'sub-' + str(pid) + '_ses-' + str(sid) + "_preproc_bold-snr-mean-416_individualspace_erode_dilate.txt")
                print("Loading %s"%filename)
                snr_mis = np.loadtxt(filename)

                # Reject bad ROIs
                rois2reject_onebased=[2, 3, 4, 6, 113, 116, 139, 203, 204, 205, 314, 319, 347]
                rois2reject_zerobased=[x-1 for x in rois2reject_onebased]
                snr_trimmed = np.delete(snr_mis, rois2reject_zerobased)
                print('The shape of the snr_trimmed array is:')
                print(snr_trimmed.shape)
                print()
                snr_all[pidind, sidind, :]=snr_trimmed   

        print('The shapeof snr_all is:')
        print(snr_all.shape)
        # np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/snr_416_erode.npy', snr_all) #Change name here!!!
        np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/snr_96_erode.npy', snr_all) #Change name here!!!

def calc_fc(run_fc, timecourse_pth, allsessid):
    # Calculating functional connectivity
    maxnsess=np.max([len(x)for x in allsessid.values()])
    allfc=np.zeros((len(allsessid), maxnsess, 387, 387))

    if run_fc:
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

                # Remove bad ROIs
                rois2reject_onebased=[2, 3, 4, 6, 113, 116, 139, 203, 204, 205, 314, 319, 347]
                rois2reject_zerobased=[x-1 for x in rois2reject_onebased]
                timecourse_trimmed = np.delete(timecourse, rois2reject_zerobased, axis=1)
                print('This is timecourse_trimmed shape:')
                print(timecourse_trimmed.shape)

                # Run filtering
                timecourse_filtered = nilearn.signal.clean(timecourse_trimmed, sessions=None, detrend=True, standardize='psc', confounds=None, low_pass=0.1, high_pass=0.01, t_r=0.392, ensure_finite=False)
                print('timecourse was filtered')

                # Calculate correlation
                fc = np.corrcoef(timecourse_filtered.T)
                print('filtered timecourse functional connectivity was correlated')

                allfc[pidind, sidind, :,:]=fc

        print('The shape of the first order correlation matrix is:')
        print(allfc.shape)
        # Saving output:
        # np.save('/dhcp/fmri_anna_graham/GKgit/finger_npy/fc_416.npy', allfc) #Change name here!!!
        np.save('/dhcp/fmri_anna_graham/GKgit/finger_npy/fc_96.npy', allfc) #Change name here!!!
    
if __name__ == '__main__':
    
    allpid, allsessid = load_arrays()
    
    snr_path = '/dhcp/scanner_positioning/rois'
    timecourse_pth = '/dhcp/fmri_anna_graham/timecourses/'
    run_snr = True ### A switch for calculating SNRHP array
    run_fc = False ### A switch for calculating FC array

    loading_SNRHP(run_snr, snr_path, allpid, allsessid)
    calc_fc(run_fc, timecourse_pth, allsessid)

    



