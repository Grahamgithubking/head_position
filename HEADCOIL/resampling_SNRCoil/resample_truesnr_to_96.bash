#!/bin/bash

#### All original code is copyright of Graham King, Trinity College Institute of Neuroscience, 12/2021.

cd /dhcp/fmri_anna_graham/GKgit/head_position/HEADCOIL/participant_files/

while IFS= read -r line
do
    SUBJ=$(echo "$line" | awk -F, '{printf "sub-"; printf $1}')
    SESS=$(echo "$line" | awk -F, '{printf "ses-"; printf $2}')
    echo $SUBJ $SESS
    
    # Extract SNR values from individual Schaefer ROIs:
    fslmeants -i  /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_preproc_bold-snr-erode1.nii.gz \
        --label=/dhcp/fmri_anna_graham/rois/${SUBJ}/${SESS}/${SUBJ}_${SESS}_schaefer_40weeks_rois.nii.gz \
        -o /dhcp/scanner_positioning/rois/${SUBJ}/${SESS}/${SUBJ}_${SESS}_snr_true.txt

done < "96twoscan.csv"
