#!/bin/bash

#### All original code is copyright of Graham King, Trinity College Institute of Neuroscience, 12/2021.

# Re-dilating the snr-mean-416-erode
# Edit sphere size in mm or box size in voxels
# Edit dilation Mean -dilM versus Mode -dilD

cd /dhcp/scanner_positioning/resliced/

# Dilated first time:
fslmaths /dhcp/scanner_positioning/resliced/preproc_bold-snr-mean-416-erode.nii.gz \
        -kernel boxv 3x3x3 -dilM \
        /dhcp/scanner_positioning/resliced/preproc_bold-snr-mean-416-erode_dilate1.nii.gz

# Dilated second time:
fslmaths /dhcp/scanner_positioning/resliced/preproc_bold-snr-mean-416-erode_dilate1.nii.gz \
        -kernel boxv 3x3x3 -dilM \
        /dhcp/scanner_positioning/resliced/preproc_bold-snr-mean-416-erode_dilate2.nii.gz


