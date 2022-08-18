from statistics import mean
import pandas as pd
import numpy as np
import glob
import os
import nilearn.signal
import matplotlib.pyplot as plt
from scipy import stats
import pingouin as pg
import seaborn as sns
from scipy.stats import bootstrap

#####################################
# arrg=np.array([[[[1],[2]],[[3],[4]]],[[[9],[10]],[[11],[12]]],[[[5],[6]],[[7],[8]]]])
# print('The shape of the array is:')
# print(arrg.shape)
# print()


# gks=np.reshape(arrg, (3*2*2,1))
# print('The gks array is:')
# print(gks)
# print()
# print('The shape of the array is:')
# print(gks.shape)


# ########################################
# x = np.array([[4,4,4],[2,2,2],[3,3,3]])
# print(x.shape)
# y=np.reshape(x, (3*3))
# print('This is y:')
# print(y)


##########################################
# x=np.array([[25.11, 30.1, np.nan, 32.02, 43.15],[14.95, 16.06, 121.25, 94.35, 29.81]])
# print('the shape of x is:')
# print(x.shape)
# print()
# y=stats.zscore(x, axis=1, nan_policy='omit')
# print(y)
# print()

################################################

# single=np.load('/dhcp/fmri_anna_graham/GKgit/finger_npy/416allsessid.npy', allow_pickle=True)
# single=dict(single)
# double=np.load('/dhcp/fmri_anna_graham/GKgit/finger_npy/96allsessid.npy', allow_pickle=True)
# double=dict(double)

# def Merge(dict1, dict2):
#     res = {**dict1, **dict2}
#     return res
     
# # Driver code
# dict1 = single
# dict2 = double
# dict3 = Merge(dict1, dict2)
# print(dict3)
# print()

# df = pd.DataFrame(single)
# dft=df.T
# dft.to_csv('/dhcp/fmri_anna_graham/GKgit/head_position/FINGER/finger_data/onesessid.csv')

########################################################


# BetaSNRcoil=np.random.rand(30,1)
# print(BetaSNRcoil)
# print()

# btpsnr=bootstrap(BetaSNRcoil.T, np.mean)
# print(btpsnr.confidence_interval)
# print(btpsnr.standard_error)

###############################################

# ofc=np.random.rand(6,2)
# ofc=np.array([[1,2],[3,4],[5,6],[7,8],[9,10]])
# print(ofc)
# print(ofc.shape)

# ofcim, pval=stats.spearmanr(ofc, axis=1)

# print(ofcim.shape)
# print()
# print(ofcim)

# plt.figure(1)
# plt.imshow(ofcim)
# plt.colorbar()
# plt.title('OFC Correlation_')
# plt.savefig('ofcim.jpg')

############################


# dhcp_root='/dhcp/dhcp_fmri_pipeline'

# # Load list of subjects
# df_subj=pd.read_csv(os.path.join(dhcp_root,'participants.tsv'), delimiter='\t')

# # Load list of sessions
# allpid=[]
# allsessid={}

# for pid in df_subj['participant_id']:
#     df_sess=pd.read_csv(os.path.join(dhcp_root, 'sub-' + pid, 'sub-' + pid + "_sessions.tsv" ), delimiter='\t')
#     if (len(df_sess))>1:
#         # More than one session
#         allpid.append(pid)
#         allsessid[pid]=[]
#         for sid in df_sess['session_id']:
#             allsessid[pid].append(sid)

# print(allpid)

#################################


mofc=np.genfromtxt('/dhcp/fmri_anna_graham/GKgit/head_position/FINGER/finger_data/48_mOFC.csv', delimiter=',')
mofc=mofc[:,1:]
mofc=np.ravel(mofc)
print(mofc)
print(mofc.shape)
np.save('/dhcp/fmri_anna_graham/GKgit/finger_npy/48_mofc.npy', mofc)
