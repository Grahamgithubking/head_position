#### All original code is copyright of Graham King, Trinity College Institute of Neuroscience.

### Plotting histograms of m-values, the slopes of the regressors

import pandas as pd
import numpy as np
import glob
import os
import nilearn.signal
import matplotlib.pyplot as plt
from scipy import stats
import pingouin as pg
import statsmodels.api as sm
import seaborn as sns

## the shape of the regressors parameters 2-D array is: [74691, 2]
## y = c + mx, where (c, m) are in parameters

parameters=np.load('/dhcp/fmri_anna_graham/GKgit/snr_npy/params_416_erode.npy') # Edit name here!
# parameters=np.load('/dhcp/fmri_anna_graham/GKgit/snr_npy/params_96_erode.npy') # Edit name here!


plt.figure(1)
sns.distplot(parameters[:,1:], kde=True) # choosing m values (aka Beta values, the slopes)
# plt.title('Histogram of the Beta values of 74,691 linear regressors \n across 416 participants')
plt.savefig('/dhcp/fmri_anna_graham/GKgit/fingerprinting/SNR/snr_figures/416_m_erode.jpg') # Edit name here!

# plt.figure(2)
# sns.distplot(parameters[:,1:], kde=True) # choosing m values (aka Beta values, the slopes)
# # plt.title('Histogram of the Beta values of 74,691 linear regressors \n across 96 sessions')
# plt.savefig('/dhcp/fmri_anna_graham/GKgit/fingerprinting/SNR/snr_figures/96_m_erode.jpg') # Edit name here!



## One Sample t-test:

ttestm=stats.ttest_1samp(parameters[:,1:], 0) # choosing m values (the slopes)

print()
print('This is the 1sample t-test statistic for m:', ttestm.statistic)
print('This is the 1sample t-test pvalue for m:', ttestm.pvalue)
print()
