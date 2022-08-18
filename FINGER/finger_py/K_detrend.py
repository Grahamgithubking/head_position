## tbc...

import pandas as pd
import numpy as np
import glob
import os
import os.path
import nilearn.signal
import nibabel as nib
import nilearn.image

    
single=np.load('/dhcp/fmri_anna_graham/GKgit/finger_npy/416allsessid.npy', allow_pickle=True)
single=dict(single)
double=np.load('/dhcp/fmri_anna_graham/GKgit/finger_npy/96allsessid.npy', allow_pickle=True)
double=dict(double)

def Merge(dict1, dict2):
    res = {**dict1, **dict2}
    return res

allsessid = Merge(single, double)
print(allsessid)
print()

# # For testing 1 session:
# allsessid={'CC00549XX22':{'157600'}}

image_loc='/dhcp/dhcp_fmri_pipeline'
save_loc='/dhcp/scanner_positioning/resliced'


def load_image(load_path, save_path, sessid):
    for pidind, (pid, sessions) in enumerate(sessid.items()):
    # for pid, sessions in allsessid.items():
        print("Working on participant %s"%pid)
        for sidind, sid in enumerate(sessions):
        # for sid in sessions:
            print("Working on session %s"%sid)
            if os.path.exists(str(save_path)+'/sub-'+str(pid)+'/ses-'+str(sid)+ '/sub-' + str(pid) + '_ses-' + str(sid) + "_bold_detrended.nii.gz"):
                print("This %s file already exists"%pidind)
            else:
                filename = os.path.join(load_path, 'sub-' + str(pid), 'ses-' + str(sid), 'func', 'sub-' + str(pid) + '_ses-' + str(sid) + "_task-rest_desc-preproc_bold.nii.gz") 
                print("Loading %s"%filename)
                image=nib.load(filename)
                ## Chose NOT to standardize as per Friedman 2006
                image_detrend=nilearn.image.clean_img(image, detrend=True, t_r=0.392, standardize=False, confounds=None, low_pass=None, high_pass=None)
                ## Saving image
                nib.save(image_detrend, str(save_path)+'/sub-'+str(pid)+'/ses-'+str(sid)+ '/sub-' + str(pid) + '_ses-' + str(sid) + "_bold_detrended.nii.gz")

load_image(image_loc, save_loc, allsessid)



# ### In terminal:
#     sbatch K_detrendbash.sh

#     squeue