#!/bin/bash

FILES=`ls -1 /dhcp/dhcp_fmri_pipeline/*/*/func/*brainmask.nii.gz | head -n 50`

mkdir -p /dhcp/scanner_positioning/resliced

for FILENAME in $FILES; do
   SUBJ=$(echo $FILENAME | cut -d"/" -f4)
   SESS=$(echo $FILENAME | cut -d"/" -f5)
   echo $SUBJ $SESS
   
   mkdir -p /dhcp/scanner_positioning/resliced/$SUBJ
   mkdir -p /dhcp/scanner_positioning/resliced/$SUBJ/$SESS

   # Now reslice this subject to template space
    flirt -in /dhcp/dhcp_fmri_pipeline/${SUBJ}/${SESS}/func/${SUBJ}_${SESS}_task-rest_desc-preproc_space-bold_brainmask.nii.gz \
            -ref /dhcp/rhodri_registration/atlases/dhcp_volume_40weeks/template_t2.nii.gz \
            -out /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_brainmask.nii.gz \
            -applyxfm

done

fslmerge -t /dhcp/scanner_positioning/resliced/mask_merged.nii.gz /dhcp/scanner_positioning/resliced/*/*/*_brainmask.nii.gz
fslmaths /dhcp/scanner_positioning/resliced/mask_merged.nii.gz -edge -thrp 50 -bin -Tmean /dhcp/scanner_positioning/resliced/mask_mean.nii.gz