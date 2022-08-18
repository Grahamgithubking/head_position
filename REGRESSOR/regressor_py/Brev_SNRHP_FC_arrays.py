## To be continued...
# Brainmasks in individual space were only eroded by a boxv of 1x1x1 only
# tbc

import pandas as pd
import numpy as np
import glob
import os
import nilearn.signal
import matplotlib.pyplot as plt
from scipy import stats
import pingouin as pg


def load_arrays(one_sessid, two_sessid):

    if one_sessid:
        allpid = np.load('/dhcp/fmri_anna_graham/GKgit/finger_npy/416allpid.npy')
        allsessid = np.load('/dhcp/fmri_anna_graham/GKgit/finger_npy/416allsessid.npy', allow_pickle=True) #I had to allow_pickle=True so as to be able to load allsessid
    if two_sessid:
        allpid = np.load('/dhcp/fmri_anna_graham/GKgit/finger_npy/48allpid.npy')
        allsessid = np.load('/dhcp/fmri_anna_graham/GKgit/finger_npy/96allsessid.npy', allow_pickle=True) #I had to allow_pickle=True so as to be able to load allsessid

    print('This is allpid:')
    print(allpid)
    print(type(allpid))
 
    allsessid = dict(allsessid)
    print('Reconverted allsessid to type dict:')
    print(allsessid)
    print(type(allsessid))
    
    return allpid, allsessid


def loading_SNRHP(run_snr, snr_path, snrcoil416, snrcoil96, snrtrue, one_sessid, two_sessid, allpid, allsessid):

    maxnsess=np.max([len(x)for x in allsessid.values()])
    snr_all=np.zeros((len(allsessid), maxnsess, 387))
    print('The original shape of snr_all zeros is:')
    print(snr_all.shape)
    print()

    if run_snr:
        for pidind, (pid, sessions) in enumerate(allsessid.items()):
        # for pid, sessions in allsessid.items():
            print("Working on participant %s"%pid)
            for sidind, sid in enumerate(sessions):
            # for sid in sessions:
                print("Working on session %s"%sid)

                if snrcoil416:
                    filename = os.path.join(snr_path, 'sub-' + str(pid), 'ses-' + str(sid), 'sub-' + str(pid) + '_ses-' + str(sid) + "_preproc_bold-snr-mean-416_individualspace_erode1.txt")
                if snrcoil96:
                    filename = os.path.join(snr_path, 'sub-' + str(pid), 'ses-' + str(sid), 'sub-' + str(pid) + '_ses-' + str(sid) + "_preproc_bold-snr-mean-96_individualspace_erode1.txt")
                if snrtrue:
                    filename = os.path.join(snr_path, 'sub-' + str(pid), 'ses-' + str(sid), 'sub-' + str(pid) + '_ses-' + str(sid) + "_snr_true.txt")    
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
        if one_sessid:
            np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/416_snr_416_erode1.npy', snr_all) #Change name here!!!
        if two_sessid:
            if snrcoil416:
                np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/96_snr_416_erode1.npy', snr_all) #Change name here!!!
            if snrcoil96:
                np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/96_snr_96_erode1.npy', snr_all) #Change name here!!!
            if snrtrue:
                np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/96_snr_true.npy', snr_all) #Change name here!!!

        

def calc_fc(run_fc, one_sessid, two_sessid, timecourse_pth, allsessid):
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
                print()

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
        if one_sessid:
            np.save('/dhcp/fmri_anna_graham/GKgit/finger_npy/fc_416.npy', allfc)
        if two_sessid:
            np.save('/dhcp/fmri_anna_graham/GKgit/finger_npy/fc_96.npy', allfc) 
        
    
if __name__ == '__main__':
    
    snr_path = '/dhcp/scanner_positioning/rois'
    timecourse_pth = '/dhcp/fmri_anna_graham/timecourses/'

    ####SWITCHES:
    run_snr = True ### A switch for calculating SNRHP array
    run_fc = False ### A switch for calculating FC array

    one_sessid=False
    two_sessid=True

    snrcoil416=False
    snrcoil96=False
    snrtrue=True


    allpid, allsessid = load_arrays(one_sessid, two_sessid)
    loading_SNRHP(run_snr, snr_path, snrcoil416, snrcoil96, snrtrue, one_sessid, two_sessid, allpid, allsessid)
    calc_fc(run_fc, one_sessid, two_sessid, timecourse_pth, allsessid)

    