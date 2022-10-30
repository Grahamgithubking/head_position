#### All original code is copyright of Graham King, Trinity College Institute of Neuroscience.
## Code to plot the mean fc for each of the 96 sessions (48 preterm, 48 term) versus PMA at scans
# Last updated 17 10 2021

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

## loading allfc:
allfc = np.load('/dhcp/fmri_anna_graham/GKgit/finger_npy/88fc.npy')
sz=allfc.shape
print('starting allfc shape is:')
print(sz)
print()
nsubj=sz[0]
nsess=sz[1]
nroi=sz[2]

allfc_reshaped=np.reshape(allfc,(nsubj, nsess, nroi*nroi))
print('The shape of allfc_reshaped is:')
print(allfc_reshaped.shape)
print()

## Get the mean fc for each session:
sessmean=allfc_reshaped.mean(axis=2)
print('The shape of sessmean is:')
print(sessmean.shape)
print()

preterm_mean=sessmean[:,0]
print('The shape of preterm_mean is:')
print(preterm_mean.shape)
term_mean=sessmean[:,1]
print('The shape of term_mean is:')
print(term_mean.shape)
acrossmean=sessmean.mean(axis=1)
print('The shape of acrossmean is:')
print(acrossmean.shape)
np.save('/dhcp/fmri_anna_graham/GKgit/finger_npy/44_mFC.npy', acrossmean)
delta_mean=term_mean - preterm_mean
print('The shape of delta_mean is:')
print(delta_mean.shape)
np.save('/dhcp/fmri_anna_graham/GKgit/finger_npy/44_deltaFC.npy', delta_mean)



## Loading in scan age data:
df_age=pd.read_csv('/dhcp/fmri_anna_graham/GKgit/head_position/FINGER/finger_data/48subjbirth_noRS.csv')
print('The imported birth dataframe is:')
print(df_age)
print()

scan1list=df_age['scan1age']
scan2list=df_age['scan2age']

plt.figure()
plt.scatter(x=(scan1list), y=(preterm_mean), c='red')
plt.scatter(x=(scan2list), y=(term_mean), c='blue')
# plt.xticks(range(0,16,2))
# plt.xlim([0, 16])
plt.title('b)', fontsize=15, fontweight="bold")
plt.xlabel('Age (PMA) at session', fontsize=13)
plt.ylabel('Fc (mean) at session', fontsize=13)
plt.legend(('Session 1', 'Session 2'), labelcolor=('red', 'blue'))
plt.savefig('/dhcp/fmri_anna_graham/GKgit/head_position/FINGER/finger_figures/fc_figures/mFC_vs_scanage.jpg')

