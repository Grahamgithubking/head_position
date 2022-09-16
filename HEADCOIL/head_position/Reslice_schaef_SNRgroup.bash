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


    # Reslice this subject's schaefer ROIs to template space:
    flirt -in /dhcp/fmri_anna_graham/rois/${SUBJ}/${SESS}/${SUBJ}_${SESS}_schaefer_40weeks_rois.nii.gz \
            -ref $REF \
            -out /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_schaefer_in_SNRgroup.nii.gz \
            -applyxfm  

            

done < "3smallheads.csv" #edit input file here
