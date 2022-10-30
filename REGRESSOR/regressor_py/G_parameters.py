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
import statsmodels.api as sm
import seaborn as sns

#####################################
# the shape of the regressors parameters 2-D array is: [74691, 2]
# y = c + mx, where (c, m) are in parameters

### For 416 participants with 1 session:
parameters=np.load('/dhcp/fmri_anna_graham/GKgit/snr_npy/74691_params_snr_416_erode1.npy') # Edit name here!
slopes=parameters[:,1:]
med=np.median(slopes)
print(f"The median of the 74691 slopes is: {med}")
print()

plt.figure(1)
sns.distplot(slopes, kde=True) # choosing m values (aka Beta values, the slopes)
# plt.title('Histogram of the Beta values of 74,691 linear regressors \n across 416 participants')
plt.xlabel('Standardized Regression Coefficient')
plt.savefig('/dhcp/fmri_anna_graham/GKgit/head_position/REGRESSOR/regressor_figures/74691_slopes_snr_416_erode1.jpg') # Edit name here!


# ## For 48 participants with 2 sessions:
# parameters=np.load('/dhcp/fmri_anna_graham/GKgit/snr_npy/74691_params_snr_416_erode1_indep88.npy') # Edit name here!
# slopes=parameters[:,1:]
# med=np.median(slopes)
# print(f"The median of the 74691 slopes is: {med}")
# print()

# plt.figure(1)
# sns.distplot(slopes, kde=True) # choosing m values (aka Beta values, the slopes)
# # plt.title('Histogram of the Beta values of 74,691 linear regressors \n across 44 independent participants')
# plt.xlabel('Standardized Regression Coefficient')
# plt.savefig('/dhcp/fmri_anna_graham/GKgit/head_position/REGRESSOR/regressor_figures/74691_slopes_snr_416_erode1_indep88.jpg') # Edit name here!

