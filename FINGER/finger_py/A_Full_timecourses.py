#### All original code is copyright of Graham King, Trinity College Institute of Neuroscience.
##First order pearson correlation of functional BOLD ROI analysis across 2300 volumes in the timecourse

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
import seaborn as sns


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
 

def pearsoncorr(timecourse_pth, allsessid, imageroot):
    # Loading resting state BOLD timecourses for the 48 participants with two sessions:
    nrois=387
    maxnsess=np.max([len(x)for x in allsessid.values()])
    allfc=np.zeros((len(allsessid), maxnsess, nrois, nrois))     # Matrix for all [48,2,387,387] first level correlations together:
    print('the shape of original allfc array is:')
    print(allfc.shape)
    

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

            # Remove bad columns
            rois2reject_onebased=[2, 3, 4, 6, 113, 116, 139, 203, 204, 205, 314, 319, 347]
            rois2reject_zerobased=[x-1 for x in rois2reject_onebased]

            timecourse_trimmed = np.delete(timecourse, rois2reject_zerobased, axis=1)

            print(timecourse_trimmed.shape)

            # Run filtering using nilearn:
            timecourse_filtered = nilearn.signal.clean(timecourse_trimmed, sessions=None, detrend=True, standardize='psc', confounds=None, low_pass=0.1, high_pass=0.01, t_r=0.392, ensure_finite=False)
            print('timecourse was filtered')

            # Calculate Pearson first order correlation across time:
            fc = np.corrcoef(timecourse_filtered.T)
            print('filtered timecourse functional connectivity was correlated')


            # Plotting RSM matrices 387x387 edges for each session
            print('here is the RSM of each session:')
            plt.imshow(fc)
            plt.colorbar()
            print('saving the png image of first level correlation:')
            imagename= str(pid) + str(sid)
            plt.savefig('imageroot' + "%s.png" %imagename)
            plt.show()

            # Matrix for all 96 first level correlations:
            allfc[pidind, sidind, :,:]=fc

    print('The shape of the 96 first order pearson correlations array is:')
    print(allfc.shape)
    print()

    # Saving output:
    np.save('/dhcp/fmri_anna_graham/GKgit/finger_npy/96fc.npy', allfc)

def session2session():
    # Session to Session analyses
    allfc = np.load('/dhcp/fmri_anna_graham/GKgit/finger_npy/96fc.npy')
    sz=allfc.shape
    print('starting allfc shape is:')
    print(sz)
    nsubj=sz[0]
    nsess=sz[1]
    nroi=sz[2]

    ##reshaping 
    allfc_reshaped=np.reshape(allfc, (nsubj * nsess, nroi * nroi))
    print('The shape of allfc_reshaped is:')
    print(allfc_reshaped.shape)

    # Histogram of First order Pearson Correlations, 96sessions 387x387 edges:
    print('the histogram of allfc_reshaped is:')
    plt.figure(1)
    sns.distplot(allfc_reshaped.ravel(), kde=True)
    plt.savefig('/dhcp/fmri_anna_graham/GKgit/fingerprinting/FINGER/finger_figures/acrosssess_figures/Hist_order1_alledges.jpg')
    plt.show()

    ## Scatter plot of mean FC across at the edge level (across all 48 participants):
    allfc_mean=np.mean(allfc, axis=0) #allfc_mean is a matrix [2,387, 387]
    allfc_mean_reshaped=np.reshape(allfc_mean,(nsess, nroi * nroi))
    print('the shape of allfc_mean_reshaped is:')
    print(allfc_mean_reshaped.shape)
    print()
    plt.figure(2)
    print('Functional Connectivity - session 1 edges vs session 2 edges \n (mean across subjects):')
    sns.kdeplot(x=allfc_mean_reshaped[0,:],y=allfc_mean_reshaped[1,:])
    plt.savefig('/dhcp/fmri_anna_graham/GKgit/fingerprinting/FINGER/finger_figures/acrosssess_figures/edges_s1vs2.jpg')
    plt.show()



if __name__ == '__main__':
    
    allpid, allsessid= get_two_session_subjects()
    print(allsessid)
    print()

    imageroot = '/dhcp/fmri_anna_graham/96timecourses/'
    timecourse_pth = '/dhcp/fmri_anna_graham/timecourses/'

    allfc = pearsoncorr(timecourse_pth, allsessid, imageroot)

    session2session()

