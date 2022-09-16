#!/bin/bash

#### All original code is copyright of Graham King, Trinity College Institute of Neuroscience, 12/2021.

#Chops out 6 motion parameters from specific volume
#Uses par_to_affine.py to translate these to FLIRT format transformation matrix (FLIRT (FMRIB's Linear Image Registration Tool))
#Reorients the Schaefer40wk template in individual functional space representing the motion at the specific volume

cd /dhcp/fmri_anna_graham/GKgit/head_position/HEADCOIL/participant_files/

while IFS= read -r line
do
    SUBJ=$(echo "$line" | awk -F, '{printf "sub-"; printf $1}')
    SESS=$(echo "$line" | awk -F, '{printf "ses-"; printf $2}')
    echo $SUBJ $SESS
    echo

    # Chopping out row 537 (mid of 0to1073) or 1763 (mid of 1227to2299) from motion.tsv file:
    # params=$(awk 'NR==537{printf "%.6f\n %.6f\n %.6f\n %.6f\n %.6f\n %.6f\n", $1, $2, $3, $4, $5, $6}' /dhcp/dhcp_fmri_pipeline/${SUBJ}/${SESS}/func/${SUBJ}_${SESS}_motion.tsv)
    params=$(awk 'NR==1763{printf "%.6f\n %.6f\n %.6f\n %.6f\n %.6f\n %.6f\n", $1, $2, $3, $4, $5, $6}' /dhcp/dhcp_fmri_pipeline/${SUBJ}/${SESS}/func/${SUBJ}_${SESS}_motion.tsv)
    echo $params


    #Calling python:
    python /dhcp/fmri_anna_graham/GKgit/head_position/HEADCOIL/par_to_affine.py /dhcp/fmri_anna_graham/rois/${SUBJ}/${SESS}/${SUBJ}_${SESS}_schaefer_40weeks_rois.nii.gz $params --invert > eventrans.mat
    

    ##Reorienting the brainmask: Choose between output 537 and 1763 here:
    # flirt -in /dhcp/fmri_anna_graham/rois/${SUBJ}/${SESS}/${SUBJ}_${SESS}_schaefer_40weeks_rois.nii.gz \
    #         -out /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_schaefer_537.nii.gz \
    #         -ref /dhcp/fmri_anna_graham/rois/${SUBJ}/${SESS}/${SUBJ}_${SESS}_schaefer_40weeks_rois.nii.gz \
    #         -interp nearestneighbour -init eventrans.mat -applyxfm
    
    flirt -in /dhcp/fmri_anna_graham/rois/${SUBJ}/${SESS}/${SUBJ}_${SESS}_schaefer_40weeks_rois.nii.gz \
            -out /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_schaefer_1763.nii.gz \
            -ref /dhcp/fmri_anna_graham/rois/${SUBJ}/${SESS}/${SUBJ}_${SESS}_schaefer_40weeks_rois.nii.gz \
            -interp nearestneighbour -init eventrans.mat -applyxfm

    echo
    echo
            
done < "96twoscan.csv" #edit input file here
