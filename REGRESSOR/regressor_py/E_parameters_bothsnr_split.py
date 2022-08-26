#### All original code is copyright of Graham King, Trinity College Institute of Neuroscience.

### Plotting histograms of m-values, the slopes of the regressors

import pandas as pd
import numpy as np
import glob
import os
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import bootstrap
import pingouin as pg
import seaborn as sns

# SWITCHES:

edgelevel=True
subjectlevel=False

orthog=True # Set this to False if doing snrfcy!

snrtrue=True
snrcoil=True

bothsess=True

bootstrapping=True

filepath='/dhcp/fmri_anna_graham/GKgit/snr_npy/'

##################### Split-Session:

if edgelevel:
    level='edgelevel'
    x='74691'
    y='298764'  # 74691 x 4directions
if subjectlevel:
    level='subjectlevel'
    x='48'
    y='192'  #48 x 4directions

if orthog:
    orthostatus='_orthog'
else:
    orthostatus='_nonorthog'

if orthog:
    p1a=np.load(str(filepath) + str(x) + '_snrboth_orthog1A.npy')
    p1b=np.load(str(filepath) + str(x) + '_snrboth_orthog1B.npy')
    p2a=np.load(str(filepath) + str(x) + '_snrboth_orthog2A.npy')
    p2b=np.load(str(filepath) + str(x) + '_snrboth_orthog2B.npy')
else:
    p1a=np.load(str(filepath) + str(x) + '_snrboth_1A.npy')
    p1b=np.load(str(filepath) + str(x) + '_snrboth_1B.npy')
    p2a=np.load(str(filepath) + str(x) + '_snrboth_2A.npy')
    p2b=np.load(str(filepath) + str(x) + '_snrboth_2B.npy')

if bothsess:
    Betasnrcoil=np.concatenate((p1a[:,1:2],p1b[:,1:2],p2a[:,1:2],p2b[:,1:2]), axis=0)
    Betasnrtrue=np.concatenate((p1a[:,2:3],p1b[:,2:3],p2a[:,2:3],p2b[:,2:3]), axis=0)
    Betafc=np.concatenate((p1a[:,3:4],p1b[:,3:4],p2a[:,3:4],p2b[:,3:4]), axis=0)
    nsess='_bothsess'


plt.figure(1)
sns.distplot(Betasnrcoil, kde=True, color='red') # choosing m values (aka Beta values, the slopes)
sns.distplot(Betasnrtrue, kde=True, color='green')
sns.distplot(Betafc, kde=True, color='blue')
plt.legend(labels=('SNR_Coil', 'SNR_True', 'FC'), labelcolor=('red', 'green', 'blue'))
plt.title(str(level) + str(orthostatus) + str(nsess) + '_all_betas_split-session')
plt.xlabel('Standardized Regression Coefficient')
ax=plt.gca()
ax.set_xlim(-2,2)
plt.grid()
plt.savefig('/dhcp/fmri_anna_graham/GKgit/head_position/REGRESSOR/regressor_figures/' + str(level) + str(orthostatus) + str(nsess) + '_all_betas_split-session.jpg')


n=np.size(Betasnrcoil, axis=0)
dof=n-1
print('The DOF is:')
print(dof)

snrcoilmean=np.mean(Betasnrcoil, axis=0)
print('The mean of the SNRCoil sample is:')
print(snrcoilmean)
snrcoilstd=np.std(Betasnrcoil, axis=0)
print('The Std of the SNRCoil sample is:')
print(snrcoilstd)
print()
snrtruemean=np.mean(Betasnrtrue, axis=0)
print('The mean of the SNRtrue sample is:')
print(snrtruemean)
snrtruestd=np.std(Betasnrtrue, axis=0)
print('the Std of the SNRtrue sample is:')
print(snrtruestd)
print()
fcmean=np.mean(Betafc, axis=0)
print('The mean of the FC sample is:')
print(fcmean)
fcstd=np.std(Betafc, axis=0)
print('The Std of the FC sample is:')
print(fcstd)
print()

## Bootstrapping:
if bootstrapping:
    btsnrcoil=bootstrap(Betasnrcoil.T, np.mean, batch=10)
    print('This is the CI of BetaSNRCoil:')
    print(btsnrcoil.confidence_interval)
    print('This is the SE of BetaSNRCoil:')
    print(btsnrcoil.standard_error)
    print()
    btsnrtrue=bootstrap(Betasnrtrue.T, np.mean, batch=10)
    print('This is the CI of BetaSNRtrue:')
    print(btsnrtrue.confidence_interval)
    print('This is the SE of BetaSNRtrue:')
    print(btsnrtrue.standard_error)
    print()
    btfc=bootstrap(Betafc.T, np.mean, batch=10)
    print('This is the CI of BetaFC:')
    print(btfc.confidence_interval)
    print('This is the SE of BetaFC:')
    print(btfc.standard_error)
    print()
