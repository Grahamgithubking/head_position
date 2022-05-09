#!/bin/bash

# Once mean SNR has been calculated across all subjects, this reslices back to indivudal space and then samples in ROIs
# Rhodri Cusack TCIN 04/2021, cusackrh@tcd.ie

FILES=`ls -1 /dhcp/dhcp_fmri_pipeline/*/*/func/*brainmask.nii.gz`  # All sessions
# FILES=`ls -1 /dhcp/dhcp_fmri_pipeline/*/*/func/*brainmask.nii.gz | head -n 5`  # Just first 5 sessions, for testing

mkdir -p /dhcp/scanner_positioning/roi

FPTH=/dhcp/dhcp_fmri_pipeline

for FILENAME in $FILES; do
   SUBJ=$(echo $FILENAME | cut -d"/" -f4)
   SESS=$(echo $FILENAME | cut -d"/" -f5)
   echo $SUBJ $SESS

    mkdir -p /dhcp/scanner_positioning/rois/${SUBJ}
    mkdir -p /dhcp/scanner_positioning/rois/${SUBJ}/${SESS}


   # Reslice mean SNR back to indiividual subject space, if it does not already exist
   if [ ! -f /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_preproc_bold-snr-mean_individualspace.nii.gz ]; then
    flirt -in /dhcp/scanner_positioning/resliced/preproc_bold-snr-mean.nii.gz \
            -ref /dhcp/dhcp_fmri_pipeline/${SUBJ}/${SESS}/func/${SUBJ}_${SESS}_task-rest_desc-preproc_space-bold_brainmask.nii.gz \
            -out /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_preproc_bold-snr-mean_individualspace.nii.gz \
            -applyxfm
   fi

    #### extract mean SNR values from rois, if it does not already exist
   if [ ! -f /dhcp/scanner_positioning/rois/${SUBJ}/${SESS}/${SUBJ}_${SESS}_preproc_bold-snr-mean_individualspace.txt ]; then
    fslmeants -i  /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_preproc_bold-snr-mean_individualspace.nii.gz \
     --label=/dhcp/fmri_anna_graham/rois/${SUBJ}/${SESS}/${SUBJ}_${SESS}_schaefer_40weeks_rois.nii.gz \
     -o /dhcp/scanner_positioning/rois/${SUBJ}/${SESS}/${SUBJ}_${SESS}_preproc_bold-snr-mean_individualspace.txt
   fi

done




