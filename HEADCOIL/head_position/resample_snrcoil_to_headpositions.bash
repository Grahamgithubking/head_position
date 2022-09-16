#!/bin/bash

#### All original code is copyright of Graham King, Trinity College Institute of Neuroscience, 12/2021.

# SNRCoil (resliced in the individual space) is sampled to give SNR values for new [schaefer40week_537 or schaefer40week_1763] head positions:
# Graham King 08/2022 TCIN

cd /dhcp/fmri_anna_graham/GKgit/head_position/HEADCOIL/participant_files/

while IFS= read -r line
do
    SUBJ=$(echo "$line" | awk -F, '{printf "sub-"; printf $1}')
    SESS=$(echo "$line" | awk -F, '{printf "ses-"; printf $2}')
    echo $SUBJ $SESS
    
    ##SNRCoil has been already sliced back to individual functional spaces previously e.g.
    ##/dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_preproc_bold-snr-mean-416_individualspace_erode1.nii.gz

    #### Extract SNRcoil values corresponding to each re-positioned schaefer40week template:
    
    fslmeants -i  /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_preproc_bold-snr-mean-416_individualspace_erode1.nii.gz \
        --label=/dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_schaefer_537.nii.gz \
        -o /dhcp/scanner_positioning/rois/${SUBJ}/${SESS}/${SUBJ}_${SESS}_SNR_HP537.txt

    fslmeants -i  /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_preproc_bold-snr-mean-416_individualspace_erode1.nii.gz \
        --label=/dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_schaefer_1763.nii.gz \
        -o /dhcp/scanner_positioning/rois/${SUBJ}/${SESS}/${SUBJ}_${SESS}_SNR_HP1763.txt

done < "96twoscan.csv"
