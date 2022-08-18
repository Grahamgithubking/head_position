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

r = np.random.rand(2,3)
print('this is r:')
print(r)
print('r has shape:')
print(r.shape)
print()
# rcut=r[:,:1]
# rcut=r[:,1:]
# print('This is rcut:')
# print(rcut)
# print()

s = np.random.rand(2,3)
print('this is s:')
print(s)
print('s has shape:')
print(s.shape)
print()

f = np.random.rand(2,3)
print('this is f:')
print(f)
print('f has shape:')
print(f.shape)
print()


for i in range(2):
    print(r[i,:])
    print()
    print(s[i,:])
    print()
    print(f[i,:])
    print()

    x=np.stack((r[i,:],s[i,:]), axis=0)
    print('This is x:')
    print(x)
    print()

    # g=f[i,:].ravel()
    # y=g.tolist()
    y=f[i,:]
    print('this is y:')
    print(y)
    print()

    X2=sm.add_constant(x.T) ### This is the key position for .T here!!
    print('This is X2 with constants:')
    print(X2)
    print()

    results=sm.OLS(y, X2).fit()
    print(results.summary())
    print()
    print()

