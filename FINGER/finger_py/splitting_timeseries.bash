#!/bin/bash

#### All original code is copyright of Graham King, Trinity College Institute of Neuroscience, 12/2021.

# Edit input file below

cd /dhcp/fmri_anna_graham/GKgit/head_position/HEADCOIL/participant_files/

while IFS= read -r line
do
    SUBJ=$(echo "$line" | awk -F, '{printf "sub-"; printf $1}')
    SESS=$(echo "$line" | awk -F, '{printf "ses-"; printf $2}')
    echo $SUBJ $SESS

    mkdir /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/detrendvolumes/
    cd /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/detrendvolumes/

    fslsplit /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_bold_detrended.nii.gz -t vol



    mkdir /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/boldvolumes/
    cd /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/boldvolumes/

    fslsplit /dhcp/dhcp_fmri_pipeline/${SUBJ}/${SESS}/func/${SUBJ}_${SESS}_task-rest_desc-preproc_bold.nii.gz -t vol
    
    

    cd /dhcp/fmri_anna_graham/GKgit/head_position/HEADCOIL/participant_files/
    VOLa=`cat volumes0to1073.csv | awk -F, '{printf $1; printf ".nii.gz "}'`
    VOLb=`cat volumes1227to2299.csv | awk -F, '{printf $1; printf ".nii.gz "}'`


    cd /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/detrendvolumes/

    fslmerge -t /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_bold_detrended_splita.nii.gz ${VOLa}
    fslmerge -t /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_bold_detrended_splitb.nii.gz ${VOLb}

    cd /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/boldvolumes/

    fslmerge -t /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_bold_splita.nii.gz ${VOLa}
    fslmerge -t /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/${SUBJ}_${SESS}_bold_splitb.nii.gz ${VOLb}



done < "96left.csv" #edit input file here as needed
