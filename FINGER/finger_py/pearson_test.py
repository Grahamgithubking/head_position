## To check how numpy correlation works
## Use https://www.socscistatistics.com/tests/pearson/default2.aspx



import numpy as np
import pandas as pd

xarr = np.array([[2,3,9,5,6],[6,1,10,13,12],[15,23,66,80,33]])
yarr = np.array([[5,6,3,2,8],[33,2,44,5,67],[20,23,65,33,22]])
print(xarr)
print(yarr)

xcor = np.corrcoef(xarr) # if 3 samples (each has 5 observations), then the resulting matrix is 3x3
xycor = np.corrcoef(xarr, yarr) # if 3 samples from x, and, 3 samples from y, then the resulting matrix is 6x6

print(f'This is the xcorr correlation: {xcor}')
print()
print(f'This is the xycorr correlation: {xycor}')


