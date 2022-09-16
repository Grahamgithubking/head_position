#### All original code is copyright of Graham King, Trinity College Institute of Neuroscience.
## Code to compare Connectome stability versus mean OFC over the two sessions.
## Last edited 20 01 2022

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
print('Ths is within:')
print(within)
print()

# Delete the participants without recorded OFC data, index positions from 0-48:
index = [2,29,30,33,35,43]
within42 = np.delete(within, index)
print('This is within42:')
print(within42)
print()

## Load Mean OFC for the 42 participants with both OFCs measured:
data = pd.read_csv('/dhcp/fmri_anna_graham/GKgit/head_position/FINGER/finger_data/48_mOFC.csv', index_col=False, header=None, names=["sujbect", "OFC"])
print('Here is the data for mean ofc:')
print(data)
print()
meanofc = data["OFC"].tolist()
print('Here is the list of Mean OFC:')
print(meanofc)
print()


plt.figure(1)
plt.scatter(x=(meanofc), y=(within42), c='grey')
# plt.xticks(range(28,36,1))
# plt.xlim([28, 36])
plt.title('a)', fontsize=15, fontweight="bold")
plt.xlabel('Mean OFC', fontsize=13)
plt.ylabel('Connectome Stability', fontsize=13)
plt.legend(['r=0.42, *p<0.005'], loc='lower right', fontsize=13)
plt.savefig('/dhcp/fmri_anna_graham/GKgit/fingerprinting/FINGER/finger_figures/acrosssess_figures/CS_vs_OFC.jpg')

pc1 = pg.corr(meanofc, within42, tail='two-sided', method='pearson')
print('The pearson correlation coefficient for Mean OFC and Within-Subject Connectome is:')
print(pc1)
print()

