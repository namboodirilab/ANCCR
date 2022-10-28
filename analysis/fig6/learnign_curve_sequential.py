from functions.load import *
from functions.plot import *
from functions.general import *
import os
import matplotlib.pyplot as plt
import numpy as np
from itertools import compress
import scipy.stats as stat
import pandas as pd
import pingouin as pg
from statsmodels.stats.anova import AnovaRM


onedrivedir = 'D:/OneDrive - UCSF'
#onedrivedir = 'D:/OneDrive - University of California, San Francisco'
directory = onedrivedir+'/Huijeong'
dathetlist = ['M2','M3','M4','M5','M6','M7','M8']
datwtlist = ['F1']
wtlist = ['F1','F2','F3','M1','M2','M3']
mouselist = ['HJ_FP_datHT_stGtACR_'+i for i in dathetlist] + ['HJ_FP_datWT_stGtACR_'+i for i in datwtlist]\
            + ['HJ_FP_WT_stGtACR_'+i for i in wtlist]
foldername = ['randomrewards','pavlovian']
daylist = []

cs1index = 15
rewardindex = 10
randomrewardindex = 7
lickindex = 5
window = [0,1500]
window_step = [-1500,0,1500]

ref = ['cs1','cs2','reward','lick','randomreward','lickcs1']
response = {}
for i in ref:
    response[i] = {}
licks = np.full((len(mouselist),1),np.nan)
for im,mousename in enumerate(mouselist):
    print(mousename)

    dfffiles, _ = findfiles(os.path.join(directory, mousename, 'pavlovian'), '.p', daylist)
    dfffiles_rr, _ = findfiles(os.path.join(directory, mousename, 'randomrewards'), '.p', daylist)
    probetestidx = [i for i, v in enumerate(dfffiles) if 'probetest' in v]
    if len(probetestidx) > 0:
        dfffiles = dfffiles[:probetestidx[0]]+dfffiles_rr[:1]

    for i in ref:
        response[i][mousename] = []

    for i,v in enumerate(dfffiles):
        matfile, _ = findfiles(os.path.dirname(v), '.mat', [])
        matfile = load_mat(matfile[0])
        dff = load_pickle(v)
        dff = dff[0]

        if i == len(dfffiles)-1:
            firstlicktimes = first_event_time_after_reference(matfile['eventlog'], lickindex, randomrewardindex, 5000)
            dopamine_rsp = [calculate_auc(dff['dff'],dff['time'],firstlicktimes,window)]
            baseline_rsp = [calculate_auc(dff['dff'],dff['time'],firstlicktimes,[-1500,0])]
            response['randomreward'][mousename].append(np.subtract(dopamine_rsp, baseline_rsp))
        else:
            firstlicktimes = first_event_time_after_reference(matfile['eventlog'], lickindex, rewardindex, 5000)
            licktimes = matfile['eventlog'][matfile['eventlog'][:, 0] == lickindex, 1]
            cuetimes = matfile['eventlog'][matfile['eventlog'][:, 0] == cs1index, 1]
            cuetimes = cuetimes[range(0, len(cuetimes), 2)]

            dopamine_rsp = [calculate_auc(dff['dff'], dff['time'], cuetimes+iw, window) for iw in window_step]+\
                           [calculate_auc(dff['dff'],dff['time'],firstlicktimes,window)]
            for ir,vr in enumerate(ref[:-3]):
                response[vr][mousename].append(np.subtract(dopamine_rsp[ir + 1], dopamine_rsp[0]))

            lick_rsp = [[np.sum(np.logical_and(licktimes>=ic+j,licktimes<ic+j+3000)) for ic in cuetimes] for j in [-3000,0]]
            lickcs1_rsp = [[np.sum(np.logical_and(licktimes >= ic + j, licktimes < ic + j + 1500)) for ic in cuetimes] for
                           j in [-1500, 0]]
            response['lick'][mousename].append(np.subtract(lick_rsp[1],lick_rsp[0]))
            response['lickcs1'][mousename].append(np.subtract(lickcs1_rsp[1],lickcs1_rsp[0]))


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

nblock = 1
clr = ['black','red']
clr_light = ['grey','pink']
cm = 1/2.54

# fig6I left
lastdaydata = {}
fig = plt.figure(figsize=(4.5*cm,3*cm))
rect = 0.35,0.2,0.53,0.68
ax = fig.add_axes(rect)
for itype in range(0, 2):
    if itype == 0:
        intype = [i for i, v in enumerate(mouselist) if np.logical_or('WT' in v, 'wt' in v)]
    else:
        intype = [i for i, v in enumerate(mouselist) if np.logical_and('HT' in v, '8' not in v)]
    cs1rsp = [flatten([movmean(v,int(np.floor(len(v)/nblock)),int(np.floor(len(v)/nblock)),0) for v in response['cs1'][mouselist[i]]]) for i in intype]
    rwrsp = [flatten([movmean(v,int(np.floor(len(v)/nblock)),int(np.floor(len(v)/nblock)),0) for v in response['reward'][mouselist[i]]]) for i in intype]

    [plt.plot(np.divide(x[:3]+x[-2:],y[0]),color=clr_light[itype],linewidth=0.35) for x,y in zip(cs1rsp,rwrsp)]
    data = [np.divide(x[:3]+x[-2:],y[0]) for x,y in zip(cs1rsp,rwrsp)]
    p, tstat, param = lnregress(np.concatenate(np.tile(np.arange(5),(len(intype),1))).reshape(-1,1), np.concatenate(data))
    lastdaydata[itype] = [x[-1] for x in data]
    plt.errorbar(range(0,5,1),[np.mean(x) for x in zip(*[np.divide(x[:3]+x[-2:], y[0]).tolist() for x, y in zip(cs1rsp, rwrsp)])],
                 [np.std(x)/np.sqrt(len(intype)) for x in zip(*[np.divide(x[:3]+x[-2:], y[0]).tolist() for x, y in zip(cs1rsp, rwrsp)])],
                 color=clr[itype], linewidth=0.5, capsize=3)
    plt.plot([-0.5,4.5],[0,0],'k:',linewidth=0.35)
stats, p = stat.ttest_ind(lastdaydata[0],lastdaydata[1])

plt.xlabel('Session')
plt.ylabel('Norm. CS1 response\n(1=reward response\n in 1st session)')
plt.ylim([-0.2,1.5])
plt.xlim([-0.5,4.5])
plt.yticks([0,0.5,1,1.5])
plt.xticks(range(0,5),['1','2','3','n-1','n'])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
fig.savefig(onedrivedir+'//figures\manuscript\dopamine_contingency//revision/fig6_new/dynamics_sessionave_'+str(nblock)+'_block.pdf',bbox_inches='tight')

# fig6I right
lastdaydata = {}
fig = plt.figure(figsize=(4.5*cm,3*cm))
rect = 0.35,0.2,0.53,0.68
ax = fig.add_axes(rect)
for itype in range(0, 2):
    if itype == 0:
        intype = [i for i, v in enumerate(mouselist) if np.logical_or('WT' in v, 'wt' in v)]
    else:
        intype = [i for i, v in enumerate(mouselist) if np.logical_and('HT' in v, '8' not in v)]
    lickrsp = [flatten([movmean(v,int(np.floor(len(v)/nblock)),int(np.floor(len(v)/nblock)),0) for v in response['lick'][mouselist[i]]]) for i in intype]
    data = [x[:3]+x[-2:] for x in lickrsp]
    lastdaydata[itype] = [x[-1] for x in data]
    p, tstat, param = lnregress(np.concatenate(np.tile(np.arange(5),(len(intype),1))).reshape(-1,1), np.concatenate(data))

    [plt.plot(x[:3]+x[-2:],color=clr_light[itype],linewidth=0.35) for x in lickrsp]
    plt.errorbar(range(0,5,1),[np.mean(x) for x in zip(*[x[:3]+x[-2:] for x in lickrsp])],
                 [np.std(x)/np.sqrt(len(intype)) for x in zip(*[x[:3]+x[-2:] for x in lickrsp])],
                 color=clr[itype], linewidth=0.5, capsize=3)
    plt.plot([-0.5,4.5],[0,0],'k:',linewidth=0.35)
stats, p = stat.ttest_ind(lastdaydata[0],lastdaydata[1])
plt.xlabel('Session')
plt.ylabel('Anticipatory licks')
plt.ylim([-2,16])
plt.xlim([-0.5,4.5])
plt.yticks([0,5,10,15])
plt.xticks(range(0,5),['1','2','3','n-1','n'])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
fig.savefig(onedrivedir+'//figures\manuscript\dopamine_contingency//revision/fig6_new/dynamics_lick_sessionave_'+str(nblock)+'_block.pdf',bbox_inches='tight')


# figS12F
fig = plt.figure(figsize=(2*cm,3*cm))
rect = 0.35,0.2,0.53,0.68
ax = fig.add_axes(rect)
dataframe = pd.DataFrame(columns = ['session','type','mouse','lick'])
for itype in range(0, 2):
    if itype == 0:
        intype = [i for i, v in enumerate(mouselist) if np.logical_or('WT' in v, 'wt' in v)]
    else:
        intype = [i for i, v in enumerate(mouselist) if np.logical_and('HT' in v, '8' not in v)]
    lickrsp = [flatten([movmean(v,int(np.floor(len(v)/nblock)),int(np.floor(len(v)/nblock)),0) for v in response['lickcs1'][mouselist[i]]]) for i in intype]

    for i,v in enumerate(lickrsp):
        dataframe = pd.concat([dataframe, pd.DataFrame({'lick': v[:1] + v[-1:], 'type': np.tile(itype,(2,1)).flatten(),
                                               'session': np.arange(0,2), 'mouse':np.tile(intype[i],(2,1)).flatten()})],ignore_index = True)

    [plt.plot(x[:1]+x[-1:],color=clr_light[itype],linewidth=0.35) for x in lickrsp]
    plt.errorbar(range(0,2,1),[np.mean(x) for x in zip(*[x[:1]+x[-1:] for x in lickrsp])],
                 [np.std(x)/np.sqrt(len(intype)) for x in zip(*[x[:1]+x[-1:] for x in lickrsp])],
                 color=clr[itype], linewidth=0.5, capsize=3)
    _, p = stat.ttest_rel([x[0] for x in lickrsp], [x[-1] for x in lickrsp])
    plt.plot([-0.5,1.5],[0,0],'k:',linewidth=0.35)
plt.xlabel('Session')
plt.ylabel('# of aniticipatory licks\n prior to CS2')
plt.ylim([-1,5])
plt.xlim([-0.5,1.5])
plt.yticks([0,2,4])
plt.xticks(range(0,2),['1','n'])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
fig.savefig(onedrivedir+'//figures\manuscript\dopamine_contingency//revision/fig6_new/anticipatory_lick_cs1_'+str(nblock)+'_block.pdf',bbox_inches='tight')

pg.mixed_anova(dv='lick', between='type', within='session', subject ='mouse',data = dataframe)
