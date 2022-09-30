#### All original code is copyright of Graham King, Trinity College Institute of Neuroscience.
## Code to ....
## Last edited 09 2022

import pandas as pd
import numpy as np
import glob
import os
import nilearn.signal
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
from scipy import stats
import pingouin as pg
import seaborn as sns


data_preterm = pd.read_csv('/dhcp/fmri_anna_graham/GKgit/head_position/HEADCOIL/head_position/48preterm.txt', sep=' ', index_col=False, header=None, names=["x", "y", "z"])
data_term = pd.read_csv('/dhcp/fmri_anna_graham/GKgit/head_position/HEADCOIL/head_position/48term.txt', sep=' ', index_col=False, header=None, names=["x", "y", "z"])

print('Here is the data_preterm for xyz coordinates:')
print(data_preterm)
print()
print('Here is the data_term for xyz coordinates:')
print(data_term)
print()

xprem = data_preterm["x"].tolist()
top_xprem = [xprem[6], xprem[36], xprem[38]]
yprem = data_preterm["y"].tolist()
top_yprem = [yprem[6], yprem[36], yprem[38]]
xterm = data_term["x"].tolist()
top_xterm = [xterm[6], xterm[36], xterm[38]]
yterm = data_term["y"].tolist()
top_yterm = [yterm[6], yterm[36], yterm[38]]
print(f"Here is the list of xprem: {xprem}")
print()
print(f"Here is the list of yprem: {yprem}")
print()
print(f"Here is the list of xterm: {xterm}")
print()
print(f"Here is the list of yterm: {yterm}")
print()
print()
print(f"The min of y is: {min(yprem)}")
print(f"The max of y is: {max(yprem)}")
print(f"The min of x is: {min(xprem)}")
print(f"The max of x is: {max(xprem)}")
print()

plt.figure(1)
# sns.kdeplot(x=(xprem), y=(yprem), color='red')
# sns.kdeplot(x=(xterm), y=(yterm), color='blue')
plt.scatter(x=(xprem), y=(yprem), c='red')
plt.scatter(x=(xterm), y=(yterm), c='blue')
plt.xticks(range(int(min(yprem)),int(max(yprem))+1,4))
plt.yticks(range(int(min(yprem)),int(max(yprem))+1,4))
plt.xlim([-15,30])
plt.ylim([-15,30])
# plt.title('Centroids of brainmasks', fontsize=15, fontweight="bold")
plt.xlabel('right - left (mm)', fontsize=13)
plt.ylabel('anterior - posterior (mm)', fontsize=13)
plt.legend(('Preterm', 'Term'), loc='upper right', fontsize=10, labelcolor=('red', 'blue'))

### Save out:
plt.savefig('/dhcp/fmri_anna_graham/GKgit/head_position/FINGER/finger_figures/fingerfigures_misc/centroidplot.jpg')

### If want to see the 3 participants who fingerprinted successfully:
# plt.scatter(x=(top_xprem), y=(top_yprem), c='yellow', marker="D") # To identfy 3 participants who fingerprinted
# plt.scatter(x=(top_xterm), y=(top_yterm), c='green', marker="D") # To identfy 3 participants who fingerprinted
# plt.savefig('/dhcp/fmri_anna_graham/GKgit/head_position/FINGER/finger_figures/fingerfigures_misc/centroidplot_tops.jpg')

