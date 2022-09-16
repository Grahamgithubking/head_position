#!/bin/bash

#### All original code is copyright of Graham King, Trinity College Institute of Neuroscience, 12/2021.

# Edit input file below

cd /dhcp/fmri_anna_graham/GKgit/head_position/HEADCOIL/participant_files/

while IFS= read -r line
do
    SUBJ=$(echo "$line" | awk -F, '{printf "sub-"; printf $1}')
    SESS=$(echo "$line" | awk -F, '{printf "ses-"; printf $2}')
    echo $SUBJ $SESS

    rm -r /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/detrendvolumes/
    rm -r /dhcp/scanner_positioning/resliced/${SUBJ}/${SESS}/boldvolumes/
  

done < "96left.csv" #edit input file here
