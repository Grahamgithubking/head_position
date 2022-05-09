#!/bin/bash

# Calculates SNR for each subject, reslices them into a common space,
# and does a "soft mean" - for each voxel, averaging across all subjects for which that voxel was inside their brainmask
# Rhodri Cusack TCIN 04/2021, cusackrh@tcd.ie

FILES=`ls -1 /dhcp/dhcp_fmri_pipeline/*/*/func/*brainmask.nii.gz`  # All sessions
# FILES=`ls -1 /dhcp/dhcp_fmri_pipeline/*/*/func/*brainmask.nii.gz | head -n 50`  # Just first 50 sessions

mkdir -p /dhcp/scanner_positioning/resliced

# Use adult T1 to define space for reslicing
REF=/usr/local/fsl/data/standard/MNI152lin_T1_1mm.nii.gz
# REF=/dhcp/rhodri_registration/atlases/dhcp_volume_40weeks/template_t2.nii.gz

for FILENAME in $FILES; do
   SUBJ=$(echo $FILENAME | cut -d"/" -f4)
   SESS=$(echo $FILENAME | cut -d"/" -f5)
   echo $SUBJ $SESS
   
   # Make directories for output
   mkdir -p /dhcp/scanner_positioning/resliced/$SUBJ
   mkdir -p /dhcp/scanner_positioning/resliced/$SUBJ/$SESS

    # Create mean of time series, if it doesn't already exist
   if [ ! -f /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_preproc_bold-mean.nii.gz ]; then
    fslmaths /dhcp/dhcp_fmri_pipeline/${SUBJ}/${SESS}/func/${SUBJ}_${SESS}_task-rest_desc-preproc_bold.nii.gz \
            -Tmean \
            /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_preproc_bold-mean.nii.gz
   fi
    # Create st. dev.  of time series, if it doesn't already exist
   if [ ! -f /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_preproc_bold-std.nii.gz ]; then
    fslmaths /dhcp/dhcp_fmri_pipeline/${SUBJ}/${SESS}/func/${SUBJ}_${SESS}_task-rest_desc-preproc_bold.nii.gz \
            -Tstd \
            /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_preproc_bold-std.nii.gz
   fi
    # Create snr, if it doesn't already exist
   if [ ! -f /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_preproc_bold-snr.nii.gz ]; then
    fslmaths /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_preproc_bold-mean.nii.gz \
            -div /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_preproc_bold-std.nii.gz \
            -mul /dhcp/dhcp_fmri_pipeline/${SUBJ}/${SESS}/func/${SUBJ}_${SESS}_task-rest_desc-preproc_space-bold_brainmask.nii.gz \
            /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_preproc_bold-snr.nii.gz
   fi

   # Reslice this subject's SNR to template space, if it does not already exist
   if [ ! -f /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_preproc_bold-snr-resliced.nii.gz ]; then
    flirt -in /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_preproc_bold-snr.nii.gz \
            -ref $REF \
            -out /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_preproc_bold-snr-resliced.nii.gz \
            -applyxfm
   fi   

   # Reslice this subject's brainmask to template space, if it does not already exist
   if [ ! -f /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_preproc_space-bold_brainmask-resliced.nii.gz ]; then
    flirt -in /dhcp/dhcp_fmri_pipeline/${SUBJ}/${SESS}/func/${SUBJ}_${SESS}_task-rest_desc-preproc_space-bold_brainmask.nii.gz \
            -ref $REF \
            -out /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_preproc_space-bold_brainmask-resliced.nii.gz \
            -applyxfm
   fi

done
# Merge the snrs and brainmasks into one file
fslmerge -t /dhcp/scanner_positioning/resliced/preproc_bold-snr-merged.nii.gz /dhcp/scanner_positioning/resliced/*/*/*_preproc_bold-snr-resliced.nii.gz   # Note this picks up everything in this directory by filter
fslmerge -t /dhcp/scanner_positioning/resliced/preproc_space-bold_brainmask-resliced.nii.gz /dhcp/scanner_positioning/resliced/*/*/*_preproc_space-bold_brainmask-resliced.nii.gz # Note this picks up everything in this directory by filter

# Calculate soft mean - for each voxel, sum of all subjects within brain mask divided by number of those subjects
fslmaths /dhcp/scanner_positioning/resliced/preproc_bold-snr-merged.nii.gz -Tmean /dhcp/scanner_positioning/resliced/preproc_bold-snr-sum.nii.gz
fslmaths /dhcp/scanner_positioning/resliced/preproc_space-bold_brainmask-resliced.nii.gz -Tmean /dhcp/scanner_positioning/resliced/preproc_space-bold_brainmask-sum.nii.gz
fslmaths /dhcp/scanner_positioning/resliced/preproc_bold-snr-sum.nii.gz -div /dhcp/scanner_positioning/resliced/preproc_space-bold_brainmask-sum.nii.gz /dhcp/scanner_positioning/resliced/preproc_bold-snr-mean.nii.gz

