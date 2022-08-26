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


# mofc=np.genfromtxt('/dhcp/fmri_anna_graham/GKgit/head_position/FINGER/finger_data/48_mOFC.csv', delimiter=',')
# mofc=mofc[:,1:]
# mofc=np.ravel(mofc)
# print(mofc)
# print(mofc.shape)
# np.save('/dhcp/fmri_anna_graham/GKgit/finger_npy/48_mofc.npy', mofc)

#####################################

# comparesess=np.array([[4,5,2,],[15,30,12],[23,20,7]])
# print(comparesess)
# print(comparesess.shape)

# #Sort columns according to how well each subject of the session0 matches each of the session1
# comparesess_sorted=np.argsort(comparesess, axis=0)
# print(comparesess_sorted)
# print()

# # #Find out the rank of the true match (in column i, where subject i ended up in the sorted ranking)
# # rankofmatch=np.where((comparesess_sorted - np.arange(3))==0)[0]
# rankofmatch=np.where(((comparesess_sorted - np.arange(3))==0))
# # rankofmatch=np.where((comparesess_sorted - np.arange(3))==0)[1]
# print(rankofmatch)
# print()

# match = np.diag(comparesess) == np.max(comparesess, axis=0)
# print('this is match:')
# print(match)
# print()


# allpid=['sub1','sub2','sub3']
# for ind, subj in enumerate(allpid):
#     print('%s\t%d'%(subj, rankofmatch[ind]))
# print()

# x=np.arange(3)
# print(x)
# print()

# g=np.array([[0,0,0],[1,2,2],[2,1,1]])
# print(g)
# k=np.array([0,1,2])
# print(k)
# print()
# p=g - k
# print(p)

##################################

# all=np.genfromtxt('/dhcp/fmri_anna_graham/GKgit/head_position/FINGER/finger_data/48subjbirth_noRS.csv', delimiter=',')
# scanone=all[1:,2:3]
# scanone=np.ravel(scanone)
# print(scanone)
# print(scanone.shape)
# np.save('/dhcp/fmri_anna_graham/GKgit/finger_npy/48_scanone.npy', scanone)

# scanint=all[1:,4:]
# scanint=np.ravel(scanint)
# print(scanint)
# print(scanint.shape)
# np.save('/dhcp/fmri_anna_graham/GKgit/finger_npy/48_scanint.npy', scanint)

#########################################################


# fc = np.arange(48)
# fc_reshaped=np.reshape(fc, (3,4,4))
# print(fc_reshaped)

# nsubj = np.arange(3)
# reg_a = [0, 2] 
# reg_b = [0, 2] 
# mat = np.empty([len(nsubj),len(reg_a), len(reg_b)])
# for subind, subj in enumerate(nsubj):
#     for inda, a in enumerate(reg_a):
#         for indb, b in enumerate(reg_b):
#             fc_temp = fc_reshaped[subj,a,b]
#             #  a * b standing in for what ever operation you wish to perform with the two regions :)
#             mat[subind,inda,indb] = fc_temp
# print(mat)

#######################################################