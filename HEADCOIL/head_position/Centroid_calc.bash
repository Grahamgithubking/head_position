#!/bin/bash

#### All original code is copyright of Graham King, Trinity College Institute of Neuroscience, 12/2021.

# Edit input file below <48preterm or 48term>
# Eidt output .txt file below <48preterm or 48term>

cd /dhcp/fmri_anna_graham/GKgit/head_position/HEADCOIL/participant_files/

while IFS= read -r line
do
    SUBJ=$(echo "$line" | awk -F, '{printf "sub-"; printf $1}')
    SESS=$(echo "$line" | awk -F, '{printf "ses-"; printf $2}')
    echo $SUBJ $SESS

    ## Note:
        # fslstats -c = output centerofgravity in mm coordinates
        # fslstats -C = output centerofgravity in voxel coordinates

    ### Tried the BrainExtractionTool BET T2w brainmask:
    # fslstats /dhcp/dhcp_anat_pipeline/${SUBJ}/${SESS}/anat/${SUBJ}_${SESS}_desc-bet_space-T2w_brainmask.nii.gz -c | tee -a /dhcp/fmri_anna_graham/GKgit/head_position/HEADCOIL/head_position/48preterm.txt

    ### Using bold_brainmask:
    fslstats /dhcp/dhcp_fmri_pipeline/${SUBJ}/${SESS}/func/${SUBJ}_${SESS}_task-rest_desc-preproc_space-bold_brainmask.nii.gz -c | tee -a /dhcp/fmri_anna_graham/GKgit/head_position/HEADCOIL/head_position/48term.txt


done < "48term.csv" #edit input file here
