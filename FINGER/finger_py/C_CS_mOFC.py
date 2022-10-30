#### All original code is copyright of Graham King, Trinity College Institute of Neuroscience.
## Code to compare Connectome stability versus mean OFC over the two sessions.
## Last edited 20 01 2022

import pingouin as pg
import pandas as pd
import numpy as np
import glob
import os
import nilearn.signal
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
from scipy import stats


## Load within-participant spearman results (main diagonal of 44x44 RSM):
within=np.load('/dhcp/fmri_anna_graham/GKgit/finger_npy/44withinfc.npy')
print('the shape of the within array is:')
print(within.shape)
print()
print('Ths is within:')
print(within)
print()

# Delete the participants without recorded OFC data, index positions from 0-43:
index_zerobased = [2,28,29,33,41]
within39 = np.delete(within, index_zerobased)
print('This is within39:')
print(within39)
print()

## Load Mean OFC for the 42 participants with both OFCs measured:
data = pd.read_csv('/dhcp/fmri_anna_graham/GKgit/head_position/FINGER/finger_data/44_mOFC.csv', index_col=False, header=None, names=["subject", "OFC"])
print('Here is the data for mean ofc:')
print(data)
print()
meanofc = data["OFC"].tolist()
print('Here is the list of Mean OFC:')
print(meanofc)
print()
meanofc = np.delete(meanofc, index_zerobased)


plt.figure(1)
plt.scatter(x=(meanofc), y=(within39), c='grey')
# plt.xticks(range(28,36,1))
# plt.xlim([28, 36])
plt.title('a)', fontsize=15, fontweight="bold")
plt.xlabel('Mean OFC', fontsize=13)
plt.ylabel('Connectome Stability', fontsize=13)
plt.legend(['r=0.38, *p<0.016'], loc='lower right', fontsize=13)
plt.savefig('/dhcp/fmri_anna_graham/GKgit/head_position/FINGER/finger_figures/acrosssess_figures/CS_vs_OFC.jpg')

pc1 = pg.corr(meanofc, within39, alternative='two-sided', method='pearson')
print('The pearson correlation coefficient for Mean OFC and Within-Subject Connectome is:')
print(pc1)
print()

#######Partial Correlation calculations:


### A.Using SNR_true values:
snr_88=np.load('/dhcp/fmri_anna_graham/GKgit/snr_npy/88_snr_true.npy')
print('the shape of the snr_88 array is:')
print(snr_88.shape)
print()
snr_mean=snr_88.mean(axis=2).mean(axis=1)
print('The shape of snrmean is:')
print(snr_mean.shape)
print()
print('All values of snrmean are:')
print(snr_mean)
print()

# ### B.Using SNR-group values:
# snr_mean=np.load('/dhcp/fmri_anna_graham/GKgit/snr_npy/88_snr_416_erode1_mean.npy')


##################################################
meansnr = np.delete(snr_mean, index_zerobased)

datacolumns = {'CS':(within39), 'mofc':(meanofc), 'msnr':(meansnr)}
df = pd.DataFrame(datacolumns)
pd.set_option('display.max_rows', df.shape[0]+1)
print(df)
print()

ppc_one = pg.partial_corr(data=df, x=('mofc'), y=('CS'), covar=['msnr'], alternative='two-sided', method='pearson')
print('The partial correlation between mOFC and CS, controlling for mSNR(true/group), is:')
print(ppc_one)
print()
