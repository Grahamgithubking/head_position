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

r = np.array([[0.55, 0.66],[0.44, 0.33]])
print('this is r:')
print(r)
print('r has shape:')
print(r.shape)
print()
s = np.array([[45,34,37],[20,23,16]])
print('this is s:')
print(s)
print('s has shape:')
print(s.shape)
print()

# z = np.log((1+r)/(1-r))/2
# print('this is z:')
# print(z)
# print()
fz = np.arctanh(r)
print('This is fz:')
print(fz)
print()

smean=np.mean(s, axis=1)
print('Ths smean is:')
print(smean)
sstd=np.std(s, axis=1)
print('The sstd is:')
print(sstd)

sz=np.zeros((2,3))
for ind in range(2):
    sz[ind,:]=(s[ind,:]-smean[ind,])/sstd[ind,]
# sz = (s-smean)/sstd
print('This is sz:')
print(sz)
print()
