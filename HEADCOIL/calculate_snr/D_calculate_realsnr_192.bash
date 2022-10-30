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


    # Create st. dev. of DETRENDED time series:
    fslmaths /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_bold_detrended_splita.nii.gz \
            -Tstd \
            /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_bold_detrended-std_splita.nii.gz
    # Create tSNR:
    fslmaths /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_preproc_bold-mean_splita.nii.gz \
            -div /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_bold_detrended-std_splita.nii.gz \
            -mul /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_task-rest_desc-preproc_space-bold_brainmask-erode1.nii.gz \
            /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_preproc_bold-snr-erode1_splita.nii.gz



    # Create st. dev. of DETRENDED time series:
    fslmaths /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_bold_detrended_splitb.nii.gz \
            -Tstd \
            /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_bold_detrended-std_splitb.nii.gz
    # Create tSNR:
    fslmaths /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_preproc_bold-mean_splitb.nii.gz \
            -div /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_bold_detrended-std_splitb.nii.gz \
            -mul /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_task-rest_desc-preproc_space-bold_brainmask-erode1.nii.gz \
            /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_preproc_bold-snr-erode1_splitb.nii.gz
            

done < "96twoscan.csv" #edit input file here
