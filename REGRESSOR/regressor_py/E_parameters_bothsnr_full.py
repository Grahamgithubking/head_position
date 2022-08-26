#### All original code is copyright of Graham King, Trinity College Institute of Neuroscience.

### Plotting histograms of m-values, the slopes of the regressors

from statistics import mean, median
import pandas as pd
import numpy as np
import glob
import os
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import bootstrap
import pingouin as pg
# import statsmodels.api as sm
import seaborn as sns

# SWITCHES:
edgelevel=True
subjectlevel=False


orthog=True # Set this to False if doing snrfcy!

snrfcy=False  # Only used when both orthog<false> and SNRCoil<True>

bootstrapping=True

filepath='/dhcp/fmri_anna_graham/GKgit/snr_npy/'

######################################
## the shape of the regressors parameters 2-D array is: [74691, 2]
## y = c + mx, where (c, m) are in parameters

# parameters=np.load('/dhcp/fmri_anna_graham/GKgit/snr_npy/74691_params_snr_416_erode1.npy') # Edit name here!

# plt.figure(1)
# sns.distplot(parameters[:,1:], kde=True) # choosing m values (aka Beta values, the slopes)
# # plt.title('Histogram of the Beta values of 74,691 linear regressors \n across 416 participants')
# plt.savefig('/dhcp/fmri_anna_graham/GKgit/head_position/REGRESSOR/regressor_figures/416_m_erode1.jpg') # Edit name here!


#######################################
## the shape of the regressors parameters 2-D array is: [48, 3]
## y = c + m1x + m2x,  where (c, m1, m2) are in parameters

## Full-Session:

if edgelevel:
    level='_edgelevel'  # 74691 x 2directions
    x='74691'
if subjectlevel:
    level='_subjectlevel'  # 48 x 2directions
    x='48' 

if orthog:
    orthostatus='_orthog'
else:
    orthostatus='_nonorthog'


if orthog:
    p1=np.load(str(filepath) + str(x) + '_params_1_snrboth_orthog.npy')
    p2=np.load(str(filepath) + str(x) + '_params_2_snrboth_orthog.npy')
else:
    p1=np.load(str(filepath) + str(x) + '_params_1_snrboth.npy')
    p2=np.load(str(filepath) + str(x) + '_params_2_snrboth.npy')

    
Betasnrcoil=np.concatenate((p1[:,1:2],p2[:,1:2]), axis=0)
Betasnrtrue=np.concatenate((p1[:,2:3],p2[:,2:3]), axis=0)
Betafc=np.concatenate((p1[:,3:4],p2[:,3:4]), axis=0) 

plt.figure(1)
sns.distplot(Betasnrcoil, kde=True, color='red') # choosing m values (aka Beta values, the slopes)
sns.distplot(Betasnrtrue, kde=True, color='green')
sns.distplot(Betafc, kde=True, color='blue')
plt.legend(('SNR_Coil', 'SNR_True', 'FC'), labelcolor=('red', 'green', 'blue'))
plt.title(str(level) + str(orthostatus) + '_all_betas_full-session') # Edit name here!
plt.xlabel('Standardized Regression Coefficient')
ax=plt.gca()
ax.set_xlim(-2,2)
plt.grid()
plt.savefig('/dhcp/fmri_anna_graham/GKgit/head_position/REGRESSOR/regressor_figures/' + str(level) + str(orthostatus) + '_all_betas.jpg') # Edit name here!


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

