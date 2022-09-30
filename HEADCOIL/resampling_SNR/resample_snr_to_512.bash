#!/bin/bash

#### All original code is copyright of Graham King, Trinity College Institute of Neuroscience, 12/2021.

cd /dhcp/fmri_anna_graham/GKgit/head_position/HEADCOIL/participant_files/

while IFS= read -r line
do
    SUBJ=$(echo "$line" | awk -F, '{printf "sub-"; printf $1}')
    SESS=$(echo "$line" | awk -F, '{printf "ses-"; printf $2}')
    echo $SUBJ $SESS
    
    # Reslice SnrCoil back to 512 sessions
    # The 512 sessions are made up of the 416 participants with one_session and 48 participants with two_session
    # Input is snr-mean-416 AKA SnrCoil

    # EDIT THIS:
        # Snr-mean-416 has been eroded by a boxv 1x1x1mm, hence "_erode1"
        # There has been no dilation of SNRCoil
    flirt -in /dhcp/scanner_positioning/resliced/preproc_bold-snr-mean-416-erode1.nii.gz -ref /dhcp/dhcp_fmri_pipeline/${SUBJ}/${SESS}/func/${SUBJ}_${SESS}_task-rest_desc-preproc_space-bold_brainmask.nii.gz -out /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_preproc_bold-snr-mean-416_individualspace_erode1.nii.gz -applyxfm

    # Extract mean SNR values from individual Schaefer ROIs:
    fslmeants -i  /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_preproc_bold-snr-mean-416_individualspace_erode1.nii.gz \
        --label=/dhcp/fmri_anna_graham/rois/${SUBJ}/${SESS}/${SUBJ}_${SESS}_schaefer_40weeks_rois.nii.gz \
        -o /dhcp/scanner_positioning/rois/${SUBJ}/${SESS}/${SUBJ}_${SESS}_preproc_bold-snr-mean-416_individualspace_erode1.txt  #Note the use of -mean-416 here!

done < "512all.csv"
