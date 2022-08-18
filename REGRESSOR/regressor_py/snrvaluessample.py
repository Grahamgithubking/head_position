import pandas as pd
import numpy as np
import glob
import os
import nilearn.signal
import matplotlib.pyplot as plt
from scipy import stats
import pingouin as pg
import statsmodels.api as sm
import seaborn as sns

#### SWITCHES:
snrmean=True
snrtrue=False

# if snrmean:
#     snr_all=np.load('/dhcp/fmri_anna_graham/GKgit/snr_npy/96_snr_416_erode1.npy')
# if snrtrue:
#     snr_all=np.load('/dhcp/fmri_anna_graham/GKgit/snr_npy/96_snr_true.npy')

# sz=snr_all.shape
# print('The shape of snr_all is:')
# print(sz)  #snr_all has shape [48, 2, 387]
# print()
# nsubj=sz[0]
# nsess=sz[1]
# nroi=sz[2]
# print('Some values of snr_all are:')
# print(snr_all[0:3,0:1,0:10])
# print()

if snrmean:
    snr_all = np.load('/dhcp/fmri_anna_graham/GKgit/snr_npy/192hp_snr.npy')
if snrtrue:
    snr_all = np.load('/dhcp/fmri_anna_graham/GKgit/snr_npy/192hp_snrtrue.npy')
    
mz=snr_all.shape
print('starting snr_all shape is:')
print(mz) #snr_all has shape [48, 2, 2, 387]
print()
nsubj=mz[0]
nsess=mz[1]
nhp=mz[2]
nroi=mz[3]
print('Some values of snr_all are:')
print(snr_all[:3,:1,:1,:10])
print()
