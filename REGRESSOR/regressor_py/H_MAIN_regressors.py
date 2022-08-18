# tbc...

import pandas as pd
import numpy as np
import glob
import os
import matplotlib.pyplot as plt
from scipy import stats
import pingouin as pg
import statsmodels.api as sm
import seaborn as sns

### SWITCHES:
full=True
split=False

# Load arrays:
if full:
    withinfc=np.load('/dhcp/fmri_anna_graham/GKgit/finger_npy/48withinfc.npy')
    snrcoil=np.load('/dhcp/fmri_anna_graham/GKgit/snr_npy/96_snr_416_erode1.npy')
    snrtrue=np.load('/dhcp/fmri_anna_graham/GKgit/snr_npy/96_snr_true.npy')
    mfwd=np.load('/dhcp/fmri_anna_graham/GKgit/finger_npy/48_mFWD.npy')
    mofc=np.load('/dhcp/fmri_anna_graham/GKgit/finger_npy/48_mofc.npy') ## contains some 'nan' values!
if split:
    ssssss

# Get means:
if full:
    snrcoil=np.mean(snrcoil, axis=(1,2))
    snrtrue=np.mean(snrtrue, axis=(1,2))
if split:
    snrcoil=np.reshape(snrcoil, (48*2,2,387))
    snrcoil=np.mean(snrcoil, axis=(2,3))
    snrtrue=np.reshape(snrtrue, (48*2,2,387))
    snrtrue=np.mean(snrtrue, axis=(2,3))
print()
print(withinfc.shape)
print(snrcoil.shape)
print(snrtrue.shape)
print(mfwd.shape)
print(mofc.shape)
print()

## Standardize:
withinfc=stats.zscore(withinfc, nan_policy='omit')
snrcoil=stats.zscore(snrcoil, nan_policy='omit')
snrtrue=stats.zscore(snrtrue, nan_policy='omit')
mfwd=stats.zscore(mfwd, nan_policy='omit')
mofc=stats.zscore(mofc, nan_policy='omit')


## Multiple Linear Regressor:
x=np.stack((snrcoil, snrtrue, mfwd, mofc), axis=0)
y=withinfc
X2 = sm.add_constant(x.T)
results = sm.OLS(y, X2, missing='drop').fit()
parameters=results.params
p= results.pvalues[1]
r2 = results.rsquared

print('Parameters array has shape:')
print(parameters.shape)
print()

print('These are the standardized beta coefficients:')
print(f"constant: {parameters[0]}")
print(f"snrcoil: {parameters[1]}")
print(f"snrtrue: {parameters[2]}")
print(f"mfwd: {parameters[3]}")
print(f"mofc: {parameters[4]}")
print()


############################

dep=withinfc
allind=np.array([snrcoil,snrtrue,mfwd,mofc])
print(allind.shape)

allnames=['SNRCoil','SNRTrue','mFWD','mOFC']
allcolours=['y','g','b','c']


def plotting(index, ind, colour, name):
    plt.figure(index)

    plt.scatter(ind, dep, color=colour)
    min_x=np.nanmin(ind)
    max_x=np.nanmax(ind)
    print(min_x)
    print(max_x)
    x=np.linspace(min_x,max_x,48)
    y=parameters[index+1] * x + parameters[0]

    plt.plot(x,y, 'r')
    plt.xlim([min_x-0.5,max_x+0.5])
    plt.title(f"Connectome Stability versus_{name}")
    plt.xlabel(f"{name}", fontsize=12)
    plt.ylabel('Connectome Stability', fontsize=12)
    plt.legend([f"Regression Coeff = {parameters[index+1]}"], loc='upper right', labelcolor='r', fontsize=10)
    plt.savefig(f"/dhcp/fmri_anna_graham/GKgit/head_position/REGRESSOR/regressor_figures/SRC_{name}.jpg")


for index in range(4):
    name=allnames[index]
    colour=allcolours[index]
    ind=allind[index,:]
    plotting(index, ind, colour, name)
    



