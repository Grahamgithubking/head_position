#!/bin/bash

#### All original code is copyright of Graham King, Trinity College Institute of Neuroscience, 12/2021.

# Edit input file below

# Use adult T1 to define space for reslicing
REF=/usr/local/fsl/data/standard/MNI152lin_T1_1mm.nii.gz

while IFS= read -r line
do
    SUBJ=$(echo "$line" | awk -F, '{printf "sub-"; printf $1}')
    SESS=$(echo "$line" | awk -F, '{printf "ses-"; printf $2}')
    echo $SUBJ $SESS

    # Create eroded brainmask of voxels 2.1mm diameter:
    fslmaths /dhcp/dhcp_fmri_pipeline/${SUBJ}/${SESS}/func/${SUBJ}_${SESS}_task-rest_desc-preproc_space-bold_brainmask.nii.gz \
            -kernel boxv 3x3x3 -ero \
            /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_task-rest_desc-preproc_space-bold_brainmask-erode.nii.gz

    # Create snr:
    fslmaths /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_preproc_bold-mean.nii.gz \
            -div /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_preproc_bold-std.nii.gz \
            -mul /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_task-rest_desc-preproc_space-bold_brainmask-erode.nii.gz \
            /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_preproc_bold-snr-erode.nii.gz

    # Reslice this subject's SNR to template space:
    flirt -in /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_preproc_bold-snr-erode.nii.gz \
            -ref $REF \
            -out /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_preproc_bold-snr-resliced-erode.nii.gz \
            -applyxfm  

   # Reslice this subject's brainmask to template space:
    flirt -in /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_task-rest_desc-preproc_space-bold_brainmask-erode.nii.gz \
            -ref $REF \
            -out /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_preproc_space-bold_brainmask-resliced-erode.nii.gz \
            -applyxfm

done < "416onescan.csv" #edit input file here
