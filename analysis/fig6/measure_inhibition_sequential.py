from functions.load import *
from functions.plot import *
import os
import matplotlib.pyplot as plt
import numpy as np
import random
import scipy.stats as stat

onedrivedir = ['D:/OneDrive - University of California, San Francisco','D:/OneDrive - UCSF']
icom = 0
directory = os.path.join(onedrivedir[icom],'Huijeong')
dathetlist = ['M2','M3','M4','M5','M6','M7']
datwtlist = ['F1']
wtlist = ['F1','F2','F3','M1','M2','M3']
mouselist = ['HJ_FP_datHT_stGtACR_'+i for i in dathetlist] + ['HJ_FP_datWT_stGtACR_'+i for i in datwtlist]\
            + ['HJ_FP_WT_stGtACR_'+i for i in wtlist]
foldername = ['randomrewards','pavlovian']
daylist = []

cs1index = 15
rewardindex = 10
bgdrewardindex = 7
rewardindex = [bgdrewardindex,rewardindex]
refindex = [bgdrewardindex,cs1index]
lickindex = 5
window = [0,500]
window_step = range(-500,1500,500)
window_baseline = [-500,-2000]

dopamine_inh = np.full((len(mouselist),len(window_step)),np.nan)
dopamine_reward = np.full((len(mouselist),2),np.nan)
for im,mousename in enumerate(mouselist):
    print(mousename)

    dfffiles, _ = findfiles(os.path.join(directory, mousename, 'pavlovian'), '.p', daylist)
    dfffiles_rr, _ = findfiles(os.path.join(directory, mousename, 'randomrewards'), '.p', daylist)
    dfffiles = dfffiles_rr[:1] + dfffiles[:1]

    for i,v in enumerate(dfffiles):
        matfile, _ = findfiles(os.path.dirname(v), '.mat', [])
        matfile = load_mat(matfile[0])
        dff = load_pickle(v)
        dff = dff[0]

        firstlicktimes = first_event_time_after_reference(matfile['eventlog'],lickindex,rewardindex[i],3000)
        reftimes = matfile['eventlog'][matfile['eventlog'][:, 0] == refindex[i], 1]
        dopamine_baseline = np.nanmean(calculate_auc(dff['dff'], dff['time'], reftimes - np.diff(window), window))
        dopamine_reward[im, i] = np.nanmean(calculate_auc(dff['dff'], dff['time'], firstlicktimes, window)) - dopamine_baseline

        if i==1:
            cstimes = matfile['eventlog'][matfile['eventlog'][:, 0] == cs1index, 1]
            cstimes = reftimes[range(1, len(cstimes), 2)]
            dopamine_inh[im,:] = [np.nanmean(calculate_auc(dff['dff'], dff['time'], cstimes+i, window))-dopamine_baseline for i in window_step]


plt.rcParams['axes.titlesize'] = 10
plt.rcParams['axes.labelsize'] = 8
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['legend.fontsize'] = 8
plt.rcParams['legend.labelspacing'] = 0.2
plt.rcParams['axes.labelpad'] = 2
plt.rcParams['axes.linewidth'] = 0.35
plt.rcParams['xtick.major.size'] = 1
plt.rcParams['xtick.major.width'] = 0.35
plt.rcParams['xtick.major.pad'] = 3
plt.rcParams['ytick.major.size'] = 1
plt.rcParams['ytick.major.width'] = 0.35
plt.rcParams['ytick.major.pad'] = 2
plt.rcParams['lines.scale_dashes'] = False
plt.rcParams['lines.dashed_pattern'] = (2, 1)
plt.rcParams['font.sans-serif'] = ['Helvetica']
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['text.color'] = 'k'


dir = 'D:\OneDrive - University of California, San Francisco//figures\manuscript\dopamine_contingency//revision//fig6_new'
cm = 1/2.54

clr = ['black','red']
clr_light = ['grey','pink']
fig = plt.figure(figsize=(4*cm,3*cm))
rect = 0.6,0.1,0.4,0.9
rect = [x*0.85 for x in rect]
data = {}
ax = fig.add_axes(rect)
for itype in range(0,2):
    if itype == 0:
        intype = [i for i, v in enumerate(mouselist) if np.logical_or('WT' in v, 'wt' in v)]
    else:
        intype = [i for i, v in enumerate(mouselist) if 'HT' in v]
    [ax.plot(dopamine_inh[im,:]/dopamine_reward[im,1],color = clr_light[itype],linewidth = 0.35) for im in intype]
    data[itype] = [dopamine_inh[im, -1] / dopamine_reward[im, 1] for im in intype]
[stats,p] = stat.ttest_ind(data[0],data[1])

for itype in range(0,2):
    if itype == 0:
        intype = [i for i, v in enumerate(mouselist) if np.logical_or('WT' in v, 'wt' in v)]
    else:
        intype = [i for i, v in enumerate(mouselist) if 'HT' in v]
    ax.errorbar(range(0,4),np.mean([dopamine_inh[im,:]/dopamine_reward[im,1] for im in intype],0),
                 np.std([dopamine_inh[im,:]/dopamine_reward[im,1] for im in intype],0)/np.sqrt(len(intype)),
                 color=clr[itype],linewidth=0.5,capsize=3)


ax.plot([-0.5,3.5],[0,0],linewidth=0.35,linestyle=':',color='k')
plt.xticks(range(0,4),labels=[-0.25,0.25,0.75,1.25],rotation=45)
plt.yticks([-1,-0.5,0,0.5])
plt.ylim([-1,0.5])
plt.xlim([-0.5,3.5])
plt.xlabel('Time from CS2 (s)')
plt.ylabel('Norm. DA response\n (1 = reward response)')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
fig.savefig(onedrivedir[icom]+'//figures\manuscript\dopamine_contingency//revision\inhibition.pdf',bbox_inches='tight')

clr = ['red','black']
clr_light = ['pink','grey']
fig = plt.figure(figsize=(4*cm,3*cm))
rect = 0.6,0.1,0.4,0.9
rect = [x*0.85 for x in rect]
ax = fig.add_axes(rect)
for itype in range(0,2):
    if itype == 0:
        intype = [i for i, v in enumerate(mouselist) if 'HT' in v]
    else:
        intype = [i for i, v in enumerate(mouselist) if np.logical_or('WT' in v, 'wt' in v)]
    ax.bar(0.5+itype,np.mean([dopamine_reward[im,1]/dopamine_reward[im,0] for im in intype]),0.6,0,color=clr_light[itype])
    [ax.scatter(random.uniform(0,0.6)+0.2+itype,dopamine_reward[im,1]/dopamine_reward[im,0],1,color=clr[itype]) for im in intype]
    ax.errorbar(0.5+itype,np.mean([dopamine_reward[im,1]/dopamine_reward[im,0] for im in intype]),
                 np.std([dopamine_reward[im,1]/dopamine_reward[im,0] for im in intype])/np.sqrt(len(mouselist)),
                color=clr[itype],capsize=3,linewidth=0.5)
    data[itype] = [dopamine_reward[im,1]/dopamine_reward[im,0] for im in intype]
[stats, p] = stat.ttest_ind(data[0], data[1])

plt.xticks([0.5,1.5], labels=['DAT-cre','WT'],rotation=45)
plt.ylabel('Norm. reward response\n(1 = random reward)')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xlim([0,2])
fig.savefig(onedrivedir[icom]+'//figures\manuscript\dopamine_contingency//revision/fig6_new//inhibition_reward.pdf',bbox_inches='tight')
