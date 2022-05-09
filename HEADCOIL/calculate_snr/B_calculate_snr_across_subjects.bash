#!/bin/bash

#### All original code is copyright of Graham King, Trinity College Institute of Neuroscience, 12/2021.

# Edit input csv file below
# Edit ouput file name below

SNR=`cat 416onescan.csv | awk -F, '{printf "sub-"; printf $1; printf "/ses-"; printf $2; printf "/sub-"; printf $1; printf "_ses-"; printf $2; printf "_preproc_bold-snr-resliced-erode.nii.gz "}'`
MASK=`cat 416onescan.csv | awk -F, '{printf "sub-"; printf $1; printf "/ses-"; printf $2; printf "/sub-"; printf $1; printf "_ses-"; printf $2; printf "_preproc_space-bold_brainmask-resliced-erode.nii.gz "}'`

cd /dhcp/scanner_positioning/resliced

# Merge the snrs and brainmasks into one file
fslmerge -t /dhcp/scanner_positioning/resliced/preproc_bold-snr-merged-erode.nii.gz ${SNR}
fslmerge -t /dhcp/scanner_positioning/resliced/preproc_space-bold_brainmask-merged-erode.nii.gz ${MASK}


# Calculate soft mean - for each voxel, sum of all subjects within brain mask divided by number of those subjects
fslmaths /dhcp/scanner_positioning/resliced/preproc_bold-snr-merged-erode.nii.gz -Tmean /dhcp/scanner_positioning/resliced/preproc_bold-snr-sum-erode.nii.gz
fslmaths /dhcp/scanner_positioning/resliced/preproc_space-bold_brainmask-merged-erode.nii.gz -Tmean /dhcp/scanner_positioning/resliced/preproc_space-bold_brainmask-sum-erode.nii.gz
fslmaths /dhcp/scanner_positioning/resliced/preproc_bold-snr-sum-erode.nii.gz -div /dhcp/scanner_positioning/resliced/preproc_space-bold_brainmask-sum-erode.nii.gz /dhcp/scanner_positioning/resliced/preproc_bold-snr-mean-416-erode.nii.gz #edit here
