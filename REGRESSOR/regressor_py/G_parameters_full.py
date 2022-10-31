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

# SWITCHES:
edgelevel=True
subjectlevel=False

snrcoil=False
snrtrue=True

orthog=True # Set this to False if doing snrfcy!

snrfcy=False  # Only used when both orthog<false> and SNRCoil<True>

bootstrapping=True

filepath='/dhcp/fmri_anna_graham/GKgit/snr_npy/'

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

if snrfcy:
    snrstatus='_SNRCoil2'
else:
    snrstatus='_SNRCoil1'

if snrcoil:
    if orthog:
        if snrfcy:
            p1mean=np.load(str(filepath) + str(x) + '_params_snrfcy1_snr_416_erode1_orthog.npy')
            p2mean=np.load(str(filepath) + str(x) + '_params_snrfcy2_snr_416_erode1_orthog.npy')
        else:
            p1mean=np.load(str(filepath) + str(x) + '_params_1_snr_416_erode1_orthog.npy')
            p2mean=np.load(str(filepath) + str(x) + '_params_2_snr_416_erode1_orthog.npy')
    else:
        if snrfcy:
            p1mean=np.load(str(filepath) + str(x) + '_params_snrfcy1_snr_416_erode1.npy')
            p2mean=np.load(str(filepath) + str(x) + '_params_snrfcy2_snr_416_erode1.npy')
        else:
            p1mean=np.load(str(filepath) + str(x) + '_params_1_snr_416_erode1.npy')
            p2mean=np.load(str(filepath) + str(x) + '_params_2_snr_416_erode1.npy')
    
    BetaSNRcoil=np.concatenate((p1mean[:,1:2],p2mean[:,1:2]), axis=0)
    BetaFC=np.concatenate((p1mean[:,2:3],p2mean[:,2:3]), axis=0) 

    plt.figure(1)
    sns.distplot(BetaSNRcoil, kde=True, color='red') # choosing m values (aka Beta values, the slopes)
    sns.distplot(BetaFC, kde=True, color='blue')
    plt.legend(('SNR_Coil', 'FC'), labelcolor=('red', 'blue'))
    plt.title(str(level) + str(orthostatus) + str(snrstatus) + '_full-session') # Edit name here!
    plt.xlabel('Standardized Regression Coefficient')
    ax=plt.gca()
    ax.set_xlim(-2,2)
    plt.grid()
    plt.savefig('/dhcp/fmri_anna_graham/GKgit/head_position/REGRESSOR/regressor_figures/' + str(level) + str(orthostatus) + str(snrstatus) + '_coil_betas.jpg') # Edit name here!

    ## One Sample t-test:
    # ttest_betasnr=stats.ttest_1samp(BetaSNRmean, popmean=0, axis=0)
    # print('This is the 1sample t-test statistic for betasnrcoil:', ttest_betasnr.statistic)
    # print('This is the 1sample t-test pvalue for betasnrcoil:', ttest_betasnr.pvalue)
    # print()

    # ttest_betafc=stats.ttest_1samp(BetaFC, popmean=0, axis=0)
    # print('This is the 1sample t-test statistic for betaFC:', ttest_betafc.statistic)
    # print('This is the 1sample t-test pvalue for betaFC:', ttest_betafc.pvalue)
    # print()

    # n=np.size(BetaSNRmean, axis=0)
    # dof=n-1
    # print('The DOF is:')
    # print(dof)

    snrmean=np.mean(BetaSNRcoil, axis=0)
    print('The mean of the SNRCoil sample is:')
    print(snrmean)
    snrstd=np.std(BetaSNRcoil, axis=0)
    print('The Std of the SNRCoil sample is:')
    print(snrstd)

    fcmean=np.mean(BetaFC, axis=0)
    print('The mean of the FC sample is:')
    print(fcmean)
    fcstd=np.std(BetaFC, axis=0)
    print('The Std of the FC sample is:')
    print(fcstd)

    ## Bootstrapping:
    if bootstrapping:
        btsnr=bootstrap(BetaSNRcoil.T, np.mean, batch=10)
        print('This is the CI of BetaSNRCoil:')
        print(btsnr.confidence_interval)
        print('This is the SE of BetaSNRCoil:')
        print(btsnr.standard_error)
        print()
        btfc=bootstrap(BetaFC.T, np.mean, batch=10)
        print('This is the CI of BetaFC:')
        print(btfc.confidence_interval)
        print('This is the SE of BetaFC:')
        print(btfc.standard_error)
        print()

if snrtrue:   
    if orthog:
        p1true=np.load(str(filepath) + str(x) + '_params_1_snr_true_orthog.npy')
        p2true=np.load(str(filepath) + str(x) + '_params_2_snr_true_orthog.npy')
    else:
        p1true=np.load(str(filepath) + str(x) + '_params_1_snr_true.npy')
        p2true=np.load(str(filepath) + str(x) + '_params_2_snr_true.npy')

    BetaSNRtrue=np.concatenate((p1true[:,1:2],p2true[:,1:2]), axis=0)
    BetaFC=np.concatenate((p1true[:,2:3],p2true[:,2:3]), axis=0)

    plt.figure(1)
    sns.distplot(BetaSNRtrue, kde=True, color='red') # choosing m values (aka Beta values, the slopes)
    sns.distplot(BetaFC, kde=True, color='blue')
    plt.legend(('SNR_True', 'FC'), labelcolor=('red', 'blue'))
    plt.title(str(level) + str(orthostatus) + '_full-session') # Edit name here!
    plt.xlabel('Standardized Regression Coefficient')
    ax=plt.gca()
    ax.set_xlim(-2,2)
    plt.grid()
    plt.savefig('/dhcp/fmri_anna_graham/GKgit/head_position/REGRESSOR/regressor_figures/' + str(level) + str(orthostatus) + '_true_betas.jpg') # Edit name here!

    # ## One Sample t-test:
    # ttest_betasnr=stats.ttest_1samp(BetaSNRtrue, popmean=0, axis=0)
    # print('This is the 1sample t-test statistic for betasnrtrue:', ttest_betasnr.statistic)
    # print('This is the 1sample t-test pvalue for betasnrtrue:', ttest_betasnr.pvalue)
    # print()

    # ttest_betafc=stats.ttest_1samp(BetaFC, popmean=0, axis=0)
    # print('This is the 1sample t-test statistic for betaFC:', ttest_betafc.statistic)
    # print('This is the 1sample t-test pvalue for betaFC:', ttest_betafc.pvalue)
    # print()

    # n=np.size(BetaSNRtrue, axis=0)
    # dof=n-1
    # print('The DOF is:')
    # print(dof)

    snrmean=np.mean(BetaSNRtrue, axis=0)
    print('The mean of the SNRtrue sample is:')
    print(snrmean)
    snrstd=np.std(BetaSNRtrue, axis=0)
    print('The Std of the SNRtrue sample is:')
    print(snrstd)

    fcmean=np.mean(BetaFC, axis=0)
    print('The mean of the FC sample is:')
    print(fcmean)
    fcstd=np.std(BetaFC, axis=0)
    print('The Std of the FC sample is:')
    print(fcstd)

    ## Bootstrapping:
    if bootstrapping:
        btsnr=bootstrap(BetaSNRtrue.T, np.mean, batch=10)
        print('This is the CI of BetaSNRtrue:')
        print(btsnr.confidence_interval)
        print('This is the SE of BetaSNRtrue:')
        print(btsnr.standard_error)
        print()
        btfc=bootstrap(BetaFC.T, np.mean, batch=10)
        print('This is the CI of BetaFC:')
        print(btfc.confidence_interval)
        print('This is the SE of BetaFC:')
        print(btfc.standard_error)
        print()
