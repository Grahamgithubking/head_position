#!/bin/bash

#### All original code is copyright of Graham King, Trinity College Institute of Neuroscience, 12/2021.

# Edit input file below

# Use adult T1 to define space for reslicing
REF=/usr/local/fsl/data/standard/MNI152lin_T1_1mm.nii.gz

cd /dhcp/fmri_anna_graham/GKgit/head_position/HEADCOIL/participant_files/

while IFS= read -r line
do
    SUBJ=$(echo "$line" | awk -F, '{printf "sub-"; printf $1}')
    SESS=$(echo "$line" | awk -F, '{printf "ses-"; printf $2}')
    echo $SUBJ $SESS

    # Create eroded brainmask of voxels 2.1mm diameter: HERE used kernel box of 1x1x1 voxel:
    fslmaths /dhcp/dhcp_fmri_pipeline/${SUBJ}/${SESS}/func/${SUBJ}_${SESS}_task-rest_desc-preproc_space-bold_brainmask.nii.gz \
            -kernel boxv 1x1x1 -ero \
            /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_task-rest_desc-preproc_space-bold_brainmask-erode1.nii.gz

    # Create st. dev. of DETRENDED time series:
    fslmaths /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_bold_detrended.nii.gz \
            -Tstd \
            /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_bold_detrended-std.nii.gz

    # Create tSNR:
    fslmaths /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_preproc_bold-mean.nii.gz \
            -div /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_bold_detrended-std.nii.gz \
            -mul /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_task-rest_desc-preproc_space-bold_brainmask-erode1.nii.gz \
            /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_preproc_bold-snr-erode1.nii.gz
            

done < "96twoscan.csv" #edit input file here
