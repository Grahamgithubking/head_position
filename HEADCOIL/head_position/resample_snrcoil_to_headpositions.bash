#!/bin/bash

#### All original code is copyright of Graham King, Trinity College Institute of Neuroscience, 12/2021.

# SNRCoil (resliced in the individual space) is sampled to give SNR values for new [schaefer40week_576 or schaefer40week_1726] head positions:
# Graham King 12/2021 TCIN

while IFS= read -r line
do
    SUBJ=$(echo "$line" | awk -F, '{printf "sub-"; printf $1}')
    SESS=$(echo "$line" | awk -F, '{printf "ses-"; printf $2}')
    echo $SUBJ $SESS
    
    ##SNRCoil has been already sliced back to individual functional spaces previously e.g.
    ##/dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_preproc_bold-snr-mean-416_individualspace_erode_dilate.nii.gz

    #### Extract SNRcoil values corresponding to each re-positioned schaefer40week template:
    
    fslmeants -i  /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_preproc_bold-snr-mean-416_individualspace_erode_dilate.nii.gz \
        --label=/dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_schaefer_576.nii.gz \
        -o /dhcp/scanner_positioning/rois/${SUBJ}/${SESS}/${SUBJ}_${SESS}_SNR_HP576.txt

    fslmeants -i  /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_preproc_bold-snr-mean-416_individualspace_erode_dilate.nii.gz \
        --label=/dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_schaefer_1726.nii.gz \
        -o /dhcp/scanner_positioning/rois/${SUBJ}/${SESS}/${SUBJ}_${SESS}_SNR_HP1726.txt

done < "96twoscan.csv"
