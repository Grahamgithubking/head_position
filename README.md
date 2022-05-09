## Last updated 24/04/2022
## All original code is copyright of Graham King, Trinity College Institute of Neuroscience. 

This repository is divided into three folders: HEADCOIL, REGRESSOR, and FINGER.

HEADCOIL folder contains code relating to:
    - participant groups (csv files)
    - Calculating the SNRCoil of the headcoil (the soft mean SNR of 416 participants with one session)
    - Reslicing SNRCoil to individual functional space
    - Resampling/Extracting SNR values using Schaefer40wk template in individual space
    - Reorienting the Schaefer40wk template in individual space using 6 motion parameters (from dHCP motion.tsv files) and FLIRT. These 6 motion parameters represent the 
    difference in head position from volume0 to volume 576 or volume 1726 of the timecourse.
    - Extracting SNRCoil values corresponding to Schaefer ROIs of these two head positions (volume 576/1726) (for 48 participants with two sessions).
    

REGRESSOR folder contains code relating to:
    - inputting predicted SNR-Headposition values into numpy arrays (original head position is volume0 of timecourse).
    - inputting of predicted SNR-Headposition values into numpy arrays for two head positions per session (vol 576, vol 1726), for all sessions of the 48 particpants with both preterm and term sessions.
    - calculating (OLS method) regressor parameters for all 74,691 edges across 416 participants with one session, and across 48 participants (independant sample) with two sessions.
    - Outputs such as histogram figures of the density of beta/slope values of the 74,691 linear regressors across both 416 participants (used to calculate SNR-Coil) and across 48 participants (independent sample) with two sessions.


FINGER folder contains code relating to:
    - calculating functional connectomes, 
    - correlating (second order correlations) across segments/sessions and participants,
    - calculating predicted functional connectivity based on SNR-Headposition values for two head positions per session (vol 576, vol 1726)
    - fingerprint identification (highest rank) procedures and significance testing (Fisher)
    - Outputs/figures of fingerprinting analyses. 




