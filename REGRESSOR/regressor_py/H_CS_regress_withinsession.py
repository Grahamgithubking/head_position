# tbc...

import pandas as pd
import numpy as np
import glob
import os
import matplotlib.pyplot as plt
from scipy import stats
import pingouin as pg
import statsmodels.api as sm
from labellines import labelLine, labelLines
import seaborn as sns



withinfcprem=np.load('/dhcp/fmri_anna_graham/GKgit/finger_npy/48withinfc_splitspreterm.npy')
withinfcterm=np.load('/dhcp/fmri_anna_graham/GKgit/finger_npy/48withinfc_splitsterm.npy')
withinfc=np.vstack((withinfcprem,withinfcterm))
withinfc=withinfc.T
print(f"This is withinfc shape: {withinfc.shape}") # This is a [48,2] array of CS for preterm and term sessions
print()

snrcoil=np.load('/dhcp/fmri_anna_graham/GKgit/snr_npy/192hp_snr.npy') # This is a [48,2,2,387] array
snrcoil=np.mean(snrcoil, axis=(2,3)) # of 4d array, this should produce a [48,2] mean per session

mfwd=np.load('/dhcp/fmri_anna_graham/GKgit/finger_npy/48_fwd_splits.npy') # This is a [48x2x2] array, the 3rd dimension is the mean_fwd splits
mfwd=np.mean(mfwd, axis=(2))

### Comparing all 96 sesssions:
withinfc=np.reshape(withinfc, (48*2))
snrcoil=np.reshape(snrcoil, (48*2))
mfwd=np.reshape(mfwd, (48*2))

### Comparing only preterm/term sessions:
# withinfc=withinfc[:,0]
# snrcoil=snrcoil[:,0]
# mfwd=mfwd[:,0]

print('These are all the shapes:')
print(withinfc.shape)
print(snrcoil.shape)
print(mfwd.shape)
# print(mfwd)

## Standardize:
withinfc=stats.zscore(withinfc, nan_policy='omit')
snrcoil=stats.zscore(snrcoil, nan_policy='omit')
mfwd=stats.zscore(mfwd, nan_policy='omit')


######### Multiple Linear Regressor:
## Not including SNRtrue here!!!!
x=np.stack((snrcoil, mfwd), axis=0)
y=withinfc
X2 = sm.add_constant(x.T)
results = sm.OLS(y, X2, missing='drop').fit()
print(results.summary())
print()
parameters=results.params
pvalues = results.pvalues
print(f"The pvalues are: {pvalues}")
r2 = results.rsquared
print(f"The r2 value is: {r2}")
print()

print('Parameters array has shape:')
print(parameters.shape)
print()

print('These are the standardized beta coefficients:')
print(f"constant: {parameters[0]}")
print(f"mean SNRcoil: {parameters[1]}")
print(f"mean FWD: {parameters[2]}")
print()


############################
#### Separate figure for each component:
    
########################################################################
### One figure for all components:

dep=withinfc
allind=np.array([snrcoil,mfwd])
print(f"This is allind shape: {allind.shape}")
print()
allnames=['Mean SNR-group','Mean FWD']
allcolours=['g','b']

plt.figure(figsize=(20,20))

def plotting(index, ind, colour, name):
    
    min_x= -1.0
    max_x= +1.0
    print(min_x)
    print(max_x)
    x=np.linspace(min_x,max_x,42)
    y=parameters[index+1] * x + parameters[0]

    plt.plot(x,y, color=colour, label=str(round(parameters[index+1],2)), linewidth='3')


for index in range(2):
    name=allnames[index]
    colour=allcolours[index]
    ind=allind[index,:]
    plotting(index, ind, colour, name)

xvals = [-0.5,-0.7,0.7,0.8,0.7]
labelLines(plt.gca().get_lines(), align=False, xvals=xvals, yoffsets=-0.03, fontsize='26')
plt.legend([f"Mean SNR-group (p<{str(round(pvalues[1],3))})", \
    f"Mean FWD (p<{str(round(pvalues[2],3))})"], loc='lower center', labelcolor=['g','b'], fontsize=24)
# plt.title(f"Relationship to Connectome Stability: \n Standardized Regression Coeffiecients", fontsize=30)
plt.xticks([])
plt.yticks([])
plt.text(-0.25, -0.2, 'Degrees of freedom = 93 \n r2 = 0.176', fontsize = 24) # DOF = 96 observations - 2 parameters -1 = 93
plt.savefig(f"/dhcp/fmri_anna_graham/GKgit/head_position/REGRESSOR/regressor_figures/SRC_ALL_split.jpg")

