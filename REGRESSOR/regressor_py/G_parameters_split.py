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

snrtrue=False  # False = SNRCoil

split=True  #  Where SNR (Coil/True) is from a 'different' split session to FCy

snrfcy=False  # Only used when both orthog<false> and SNRCoil<True>
            # Here, SNR-Coil and FCy are from the same split

bothsess=True
preterm=False
term=False

bootstrapping=False

filepath='/dhcp/fmri_anna_graham/GKgit/snr_npy/'

## Split-Session:

if edgelevel:
    level='edgelevel'
    x='74691'
    y='298764'  # 74691 x 4directions
if subjectlevel:
    level='subjectlevel'
    x='48'
    y='192'  #48 x 4directions

if snrtrue:
    if orthog:
        p1a=np.load(str(filepath) + str(x) + '_true_orthog1A.npy')
        p1b=np.load(str(filepath) + str(x) + '_true_orthog1B.npy')
        p2a=np.load(str(filepath) + str(x) + '_true_orthog2A.npy')
        p2b=np.load(str(filepath) + str(x) + '_true_orthog2B.npy')
    else:
        p1a=np.load(str(filepath) + str(x) + '_true_1A.npy')
        p1b=np.load(str(filepath) + str(x) + '_true_1B.npy')
        p2a=np.load(str(filepath) + str(x) + '_true_2A.npy')
        p2b=np.load(str(filepath) + str(x) + '_true_2B.npy')
else:
    if split:
        if orthog:
            p1a=np.load(str(filepath) + str(x) + '_params_orthog1A.npy')
            p1b=np.load(str(filepath) + str(x) + '_params_orthog1B.npy')
            p2a=np.load(str(filepath) + str(x) + '_params_orthog2A.npy')
            p2b=np.load(str(filepath) + str(x) + '_params_orthog2B.npy')
        else:
            p1a=np.load(str(filepath) + str(x) + '_params_1A.npy')
            p1b=np.load(str(filepath) + str(x) + '_params_1B.npy')
            p2a=np.load(str(filepath) + str(x) + '_params_2A.npy')
            p2b=np.load(str(filepath) + str(x) + '_params_2B.npy')
    if snrfcy:
        if orthog:
            p1a=np.load(str(filepath) + str(x) + '_params_orthog1snrfcya.npy')
            p1b=np.load(str(filepath) + str(x) + '_params_orthog1snrfcyb.npy')
            p2a=np.load(str(filepath) + str(x) + '_params_orthog2snrfcya.npy')
            p2b=np.load(str(filepath) + str(x) + '_params_orthog2snrfcyb.npy')
        else:
            p1a=np.load(str(filepath) + str(x) + '_params_1snrfcya.npy')
            p1b=np.load(str(filepath) + str(x) + '_params_1snrfcyb.npy')
            p2a=np.load(str(filepath) + str(x) + '_params_2snrfcya.npy')
            p2b=np.load(str(filepath) + str(x) + '_params_2snrfcyb.npy')


if bothsess:
    BetaSNR=np.concatenate((p1a[:,1:2],p1b[:,1:2],p2a[:,1:2],p2b[:,1:2]), axis=0)
    BetaFC=np.concatenate((p1a[:,2:3],p1b[:,2:3],p2a[:,2:3],p2b[:,2:3]), axis=0)
    nsess='_bothsess'
if term:
    # Term session only:
    BetaSNR=np.concatenate((p2a[:,1:2],p2b[:,1:2]), axis=0)
    BetaFC=np.concatenate((p2a[:,2:3],p2b[:,2:3]), axis=0)
    nsess='_term'
if preterm:
    # Preterm session only:
    BetaSNR=np.concatenate((p1a[:,1:2],p1b[:,1:2]), axis=0)
    BetaFC=np.concatenate((p1a[:,2:3],p1b[:,2:3]), axis=0)
    nsess='_preterm'


if snrtrue:
    snrstatus='_SNRTrue'
else:
    snrstatus='_SNRCoil'

if orthog:
    orthostatus='_orthog'
else:
    orthostatus='_nonorthog'

if snrfcy:
    snrnum='_SNR2'
else:
    snrnum='_SNR1'


plt.figure(1)
sns.distplot(BetaSNR, kde=True, color='red') # choosing m values (aka Beta values, the slopes)
sns.distplot(BetaFC, kde=True, color='blue')
plt.legend(labels=('SNR', 'FC'), labelcolor=('red', 'blue'))
plt.title(str(level) + str(snrstatus) + str(snrnum) + str(orthostatus) + str(nsess) + '_split-session')
plt.xlabel('Standardized Regression Coefficient')
ax=plt.gca()
ax.set_xlim(-2,2)
plt.grid()
plt.savefig('/dhcp/fmri_anna_graham/GKgit/head_position/REGRESSOR/regressor_figures/' + str(level) + str(snrstatus) + str(snrnum) + str(orthostatus) +str(nsess) + '.jpg')

## One Sample t-test:
# ttest_betasnr=stats.ttest_1samp(BetaSNR, popmean=0, axis=0)
# print('This is the 1sample t-test statistic for betasnr:', ttest_betasnr.statistic)
# print('This is the 1sample t-test pvalue for betasnr:', ttest_betasnr.pvalue)
# print()
# ttest_betafc=stats.ttest_1samp(BetaFC, popmean=0, axis=0)
# print('This is the 1sample t-test statistic for betafc:', ttest_betafc.statistic)
# print('This is the 1sample t-test pvalue for betafc:', ttest_betafc.pvalue)
# print()


n=np.size(BetaSNR, axis=0)
dof=n-1
print('The DOF is:')
print(dof)

snrmean=np.mean(BetaSNR, axis=0)
print('The mean of the SNR sample is:')
print(snrmean)
snrstd=np.std(BetaSNR, axis=0)
print('The Std of the SNR sample is:')
print(snrstd)

fcmean=np.mean(BetaFC, axis=0)
print('The mean of the FC sample is:')
print(fcmean)
fcstd=np.std(BetaFC, axis=0)
print('The Std of the FC sample is:')
print(fcstd)

## Bootstrapping:
if bootstrapping:
    btsnr=bootstrap(BetaSNR.T, np.mean, batch=10)
    print('This is the CI of BetaSNR:')
    print(btsnr.confidence_interval)
    print('This is the SE of BetaSNR:')
    print(btsnr.standard_error)
    print()
    btfc=bootstrap(BetaFC.T, np.mean, batch=10)
    print('This is the CI of BetaFC:')
    print(btfc.confidence_interval)
    print('This is the SE of BetaFC:')
    print(btfc.standard_error)
    print()
