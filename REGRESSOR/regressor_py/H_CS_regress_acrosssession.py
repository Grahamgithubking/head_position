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

### SWITCHES:
full=True
split=False

# Load arrays:
if full:
    withinfc=np.load('/dhcp/fmri_anna_graham/GKgit/finger_npy/44withinfc.npy')
    snrcoil=np.load('/dhcp/fmri_anna_graham/GKgit/snr_npy/88_snr_416_erode1.npy')
    snrtrue=np.load('/dhcp/fmri_anna_graham/GKgit/snr_npy/88_snr_true.npy')
    mfwd=np.load('/dhcp/fmri_anna_graham/GKgit/finger_npy/44_mFWD.npy')
    scanint=np.load('/dhcp/fmri_anna_graham/GKgit/finger_npy/44_scanint.npy')
    scanone=np.load('/dhcp/fmri_anna_graham/GKgit/finger_npy/44_scanone.npy')
    mofc=np.load('/dhcp/fmri_anna_graham/GKgit/finger_npy/44_mofc.npy') ## contains some 'nan' values!
if split:
    ssssss

# Get means where needed:
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
print(scanint.shape)
print(scanone.shape)
print(mofc.shape)
print()

## Standardize:
withinfc=stats.zscore(withinfc, nan_policy='omit')
snrcoil=stats.zscore(snrcoil, nan_policy='omit')
snrtrue=stats.zscore(snrtrue, nan_policy='omit')
mfwd=stats.zscore(mfwd, nan_policy='omit')
scanint=stats.zscore(scanint, nan_policy='omit')
scanone=stats.zscore(scanone, nan_policy='omit')
mofc=stats.zscore(mofc, nan_policy='omit')


######### Multiple Linear Regressor:
## Not including SNRtrue here!!!!
x=np.stack((snrcoil, mfwd, scanint, scanone, mofc), axis=0)
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
print(f"scanint: {parameters[3]}")
print(f"scanone: {parameters[4]}")
print(f"mean OFC: {parameters[5]}")
print()


############################
#### Separate figure for each component:

# dep=withinfc

# top_dep=[withinfc[6], withinfc[36], withinfc[38]]

# allind=np.array([snrcoil,snrtrue,mfwd,scanint,scanone,mofc])
# print(f"This is allind shape: {allind.shape}")
# print()

# allnames=['Mean SNRCoil','Mean SNRTrue','Mean FWD','Scan interval(wk)','Scan_1 age(PMA)','Mean OFC']
# allcolours=['y','g','b','c','m','k']


# def plotting(index, ind, colour, name):
#     plt.figure(index)
#     top_ind=[ind[6],ind[36],ind[38]]

#     plt.scatter(ind, dep, marker="*", color=colour)
#     plt.scatter(top_ind, top_dep, marker="D", color=colour)
    
#     min_x=np.nanmin(ind)
#     max_x=np.nanmax(ind)
#     print(min_x)
#     print(max_x)
#     x=np.linspace(min_x,max_x,48)
#     y=parameters[index+1] * x + parameters[0]

#     plt.plot(x,y, 'r')
#     plt.xlim([min_x-0.5,max_x+0.5])
#     plt.title(f"Connectome Stability versus_{name}")
#     plt.xlabel(f"{name}", fontsize=12)
#     plt.ylabel('Connectome Stability', fontsize=12)
#     plt.legend([f"Coefficient = {parameters[index+1]}"], loc='lower right', labelcolor='r', fontsize=10)
#     plt.savefig(f"/dhcp/fmri_anna_graham/GKgit/head_position/REGRESSOR/regressor_figures/SRC_{name}.jpg")


# for index in range(6):
#     name=allnames[index]
#     colour=allcolours[index]
#     ind=allind[index,:]
#     plotting(index, ind, colour, name)
    
########################################################################
### One figure for all components:

dep=withinfc
allind=np.array([snrcoil,mfwd,scanint,scanone,mofc])
print(f"This is allind shape: {allind.shape}")
print()
allnames=['Mean SNR-group','Mean FWD','Scan interval(wk)','Scan_1 age(PMA)','Mean Head-size']
allcolours=['g','b','c','m','r']

plt.figure(figsize=(20,20))

def plotting(index, ind, colour, name):
    
    min_x= -1.0
    max_x= +1.0
    print(min_x)
    print(max_x)
    x=np.linspace(min_x,max_x,42)
    y=parameters[index+1] * x + parameters[0]

    plt.plot(x,y, color=colour, label=str(round(parameters[index+1],2)), linewidth='3')


for index in range(5):
    name=allnames[index]
    colour=allcolours[index]
    ind=allind[index,:]
    plotting(index, ind, colour, name)

xvals = [-0.5,-0.7,0.7,0.8,0.7]
labelLines(plt.gca().get_lines(), align=False, xvals=xvals, yoffsets=-0.03, fontsize='26')
plt.legend([f"Mean SNR-group (p<{str(round(pvalues[1],3))})", \
    f"Mean FWD (p<{str(round(pvalues[2],3))})", f"Scan interval (p<{str(round(pvalues[3],3))})", \
        f"Scan_1 PMA (p<{str(round(pvalues[4],3))})", f"Mean Head-size (p<{str(round(pvalues[5],3))})"], \
            loc='lower center', labelcolor=['g','b','c','m','r'], fontsize=24)
# plt.title(f"Relationship to Connectome Stability: \n Standardized Regression Coeffiecients", fontsize=30)
plt.xticks([])
plt.yticks([])
plt.text(-0.25, -0.2, 'Degrees of freedom = 33 \n r2 = 0.49', fontsize = 24) # DOF = 39 observations - 5 parameters -1 = 33
plt.savefig(f"/dhcp/fmri_anna_graham/GKgit/head_position/REGRESSOR/regressor_figures/SRC_ALL.jpg")

