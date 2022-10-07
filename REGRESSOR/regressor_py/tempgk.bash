#!/bin/bash

#### All original code is copyright of Graham King, Trinity College Institute of Neuroscience, 12/2021.


cd /dhcp/fmri_anna_graham/GKgit/head_position/HEADCOIL/participant_files

while IFS= read -r line
do
    SUBJ=$(echo "$line" | awk -F, '{printf "sub-"; printf $1}')
    SESS=$(echo "$line" | awk -F, '{printf "ses-"; printf $2}')
    echo $SUBJ $SESS

    # cp /dhcp/dhcp_fmri_pipeline/${SUBJ}/${SESS}/func/${SUBJ}_${SESS}_task-rest_desc-preproc_space-bold_brainmask.nii.gz /dhcp/fmri_anna_graham/GKgit/head_position/tempgkfolder/
    cp /dhcp/fmri_anna_graham/rois/${SUBJ}/${SESS}/${SUBJ}_${SESS}_schaefer_40weeks_rois.nii.gz /dhcp/fmri_anna_graham/GKgit/head_position/tempgkfolder/
    
done < "96twoscan.csv" #edit input file here

