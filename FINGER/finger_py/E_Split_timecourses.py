#### All original code is copyright of Graham King, Trinity College Institute of Neuroscience.
##First order pearson correlation of functional BOLD ROI analysis across Split Segments timecourse


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


def pearsoncorr(split_pth, allsessid, imageroot):
    # Loading 192 split segments timecourse numpy files
    nrois=387
    nsplit=2
    nsess=2
    npart=48
    allfc=np.zeros((npart, nsess, nsplit, nrois, nrois))
    print('the shape of original allfc array is:')
    print(allfc.shape)

    splits = ['_A', '_B']

    for pidind, (pid, sessions) in enumerate(allsessid.items()):
    # for pid, sessions in allsessid.items():
        print("Working on participant %s"%pid)
        for sidind, sid in enumerate(sessions):
        # for sid in sessions:
            print("Working on session %s"%sid)
            for splitind, split in enumerate(splits):
            # for split in splits:
                print("Working on split %s"%split)

                # Load up timecourse
                timecourse = np.load(split_pth + str(pid) + str(sid) + str(split) + ".npy")
                print('Loading split segment timecourse')
                print(timecourse.shape)

                # Remove bad columns
                rois2reject_onebased=[2, 3, 4, 6, 113, 116, 139, 203, 204, 205, 314, 319, 347]
                rois2reject_zerobased=[x-1 for x in rois2reject_onebased]

                timecourse_trimmed = np.delete(timecourse, rois2reject_zerobased, axis=1)
                print(timecourse_trimmed.shape)

                # Run filtering using nilearn:
                timecourse_filtered = nilearn.signal.clean(timecourse_trimmed, sessions=None, detrend=True, standardize='psc', confounds=None, low_pass=0.1, high_pass=0.01, t_r=0.392, ensure_finite=False)
                print('timecourse was filtered')

                # Calculate Pearson first order correlation across split segment:
                fc = np.corrcoef(timecourse_filtered.T)
                print('The filtered split segment timecourse functional connectivity was correlated')
           
                # Plotting RSM matrices 387x387 edges for each split segment:
                print('here is the RSM of each split:')
                plt.imshow(fc)
                plt.colorbar()
                print('saving the png image of first level correlation:')
                imagename= str(pid) + str(sid) + str(split)
                plt.savefig(imageroot + "%s.png" %imagename)
                plt.show()

                # Matrix for all 192 first level correlations:
                allfc[pidind, sidind, splitind, :,:]=fc

    print('The shape of the 192 first level correlations array is:')
    print(allfc.shape)
    print()

    # Saving output:
    np.save('/dhcp/fmri_anna_graham/GKgit/finger_npy/192fc.npy', allfc)



if __name__ == '__main__':
    
    allpid, allsessid= get_two_session_subjects()
    print(allsessid)

    split_pth = '/dhcp/fmri_anna_graham/timecourses_split/'
    imageroot = '/dhcp/fmri_anna_graham/gk_fc_192/'
    
    pearsoncorr(split_pth, allsessid, imageroot)

    
