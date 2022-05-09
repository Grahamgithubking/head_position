#### All original code is copyright of Graham King, Trinity College Institute of Neuroscience.
## Code to compare Connectome stability versus mean Predicted SNR over the two sessions.
## Last edited 16 10 2021

import pandas as pd
import numpy as np
import glob
import os
import nilearn.signal
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
from scipy import stats
import pingouin as pg


## Load within-participant spearman results (main diagonal of 48x48 RSM):
within=np.load('/dhcp/fmri_anna_graham/GKgit/finger_npy/48withinfc.npy')
print('the shape of the within array is:')
print(within.shape)
print()

## Load Predicted SNR for the 96 sessions:
snr_96=np.load('/dhcp/fmri_anna_graham/GKgit/snr_npy/snr_96_erode.npy')
print('the shape of the snr_96 array is:')
print(snr_96.shape)
print()

snr_mean=snr_96.mean(axis=2).mean(axis=1)
print('The shape of snrmean is:')
print(snr_mean.shape)
print()
print('All values of snrmean are:')
print(snr_mean)
print()
np.save('/dhcp/fmri_anna_graham/GKgit/snr_npy/snrmismean.npy', snr_mean)

plt.figure(1)
plt.scatter(x=(snr_mean), y=(within), c='grey')
# plt.xticks(range(0,16,2))
# plt.xlim([0, 16])
plt.title('b)', fontsize=15, fontweight="bold")
plt.xlabel('Mean Predicted SNR', fontsize=13)
plt.ylabel('Connectome Stability', fontsize=13)
plt.legend(['r=0.5, *p<0.0003'], loc='lower right', fontsize=13)
plt.savefig('/dhcp/fmri_anna_graham/GKgit/fingerprinting/FINGER/finger_figures/acrosssess_figures/CS_vs_mSNR.jpg')

pc1 = pg.corr(snr_mean, within, tail='two-sided', method='pearson')
print('The pearson correlation coefficient for Snrmismean and Within-Subject Connectome is:')
print(pc1)
print()

