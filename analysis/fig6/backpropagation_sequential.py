from functions.load import *
from functions.plot import *
from functions.general import *
import os
import matplotlib.pyplot as plt
import numpy as np
import random

directory = 'D:/OneDrive - UCSF/Huijeong'
#directory = 'D:\OneDrive - University of California, San Francisco\Huijeong'
dathetlist = ['M2','M3','M4','M5','M6','M7']
datwtlist = ['F1']
wtlist = ['F1','F2','F3','M1','M2','M3']
mouselist = ['HJ_FP_datWT_stGtACR_'+i for i in datwtlist] + ['HJ_FP_WT_stGtACR_'+i for i in wtlist]
foldername = ['randomrewards','pavlovian']
daylist = []

cs1index = 15
rewardindex = 10
lickindex = 5
window = [0, 3000]
window_bl = [-2000, 0]
binsize = 20

peakidx = {}
auc = {}
for im, mousename in enumerate(mouselist):
    print(mousename)

    dfffiles, _ = findfiles(os.path.join(directory, mousename, 'pavlovian'), '.p', daylist)
    probetestidx = [i for i, v in enumerate(dfffiles) if 'probetest' in v]
    if len(probetestidx) > 0:
        dfffiles = dfffiles[:probetestidx[0]]

    peakidx[mousename] = []
    auc[mousename] = {}
    auc[mousename]['early'] = []
    auc[mousename]['late'] = []
    for i, v in enumerate(dfffiles):
        matfile, _ = findfiles(os.path.dirname(v), '.mat', [])
        matfile = load_mat(matfile[0])
        dff = load_pickle(v)
        dff = dff[0]

        signal, _, _, time = align_signal_to_reference(dff['dff'], dff['time'], matfile['eventlog'], cs1index, [-2000, 3000], 0, 0)
        cs1time = matfile['eventlog'][matfile['eventlog'][:,0]==cs1index,1]
        early_temp = calculate_auc(dff['dff'], dff['time']/1000, cs1time[range(0,len(cs1time),2)]/1000, [0,1])
        late_temp = calculate_auc(dff['dff'],dff['time']/1000,cs1time[range(1,len(cs1time),2)]/1000,[0,1])
        if i==0:
            firstlicktimes = first_event_time_after_reference(matfile['eventlog'], lickindex, rewardindex, 5000)
            auc[mousename]['reward'] = calculate_auc(dff['dff'], dff['time'] / 1000, firstlicktimes/ 1000, [0, 1])
        auc[mousename]['early'] = auc[mousename]['early'] + early_temp
        auc[mousename]['late'] = auc[mousename]['late'] + late_temp

        mov_signal = movmean(signal[:, np.logical_and(time >= 0, time <= 3000)], int(binsize / np.mean(np.diff(dff['time']))), 1, 1)
        mov_baseline = movmean(signal[:, np.logical_and(time >= -2000, time < 0)], int(binsize / np.mean(np.diff(dff['time']))), 1, 1)
        mov_time = movmean(time[np.logical_and(time >= 0, time <= 3000)], int(binsize / np.mean(np.diff(dff['time']))), 1, 0)

        time_signal = time[np.logical_and(time >=0, time <= 3000)]
        peakidx_temp, _ = peaksearch(mov_signal, [2*np.std(i) for i in mov_baseline], 'max')

        peakidx[mousename] = peakidx[mousename] + peakidx_temp

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


ntrials = 200
#dir = 'D:\OneDrive - University of California, San Francisco//figures\manuscript\dopamine_contingency//revision//fig6_new'
dir = 'D:\OneDrive - UCSF//figures\manuscript\dopamine_contingency//revision//fig6_new'

cm = 1/2.54
fig = plt.figure(figsize=(4*cm, 3*cm))
rect = 0.2,0.2,0.68,0.68
ax = fig.add_axes(rect)
[ax.plot(range(0,ntrials),np.divide(auc[x]['early'][:ntrials],np.nanmean(auc[x]['early'][ntrials-50:ntrials])),'grey',linewidth=0.35) for x in mouselist]
[ax.plot(range(0,ntrials),np.divide(auc[x]['late'][:ntrials],np.nanmean(auc[x]['early'][ntrials-50:ntrials])),'lightblue',linewidth=0.35) for x in mouselist]
ax.plot(range(0,ntrials),np.nanmean([np.divide(auc[x]['early'][:ntrials],np.nanmean(auc[x]['early'][ntrials-50:ntrials])) for x in mouselist],axis=0),'black',linewidth=1)
ax.plot(range(0,ntrials),np.nanmean([np.divide(auc[x]['late'][:ntrials],np.nanmean(auc[x]['early'][ntrials-50:ntrials])) for x in mouselist],axis=0),'blue',linewidth=1)
plt.yticks(range(-2,4,1))
plt.xticks(range(0,ntrials+1,50),rotation=45)
plt.xlabel('Trial')
plt.ylabel('Norm. DA response')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xlim([0,ntrials])
plt.ylim([-1.5,2.5])
fig.savefig(dir+'//backpropagation_2cues_timecourse.pdf',bbox_inches='tight')

fig = plt.figure(figsize=(2.25*cm, 3*cm))
rect = 0.6,0.1,0.4,0.9
rect = [x*0.85 for x in rect]
ax = fig.add_axes(rect)
data = [np.divide(np.nanmean(np.subtract(auc[x]['late'][:50],auc[x]['early'][:50])),np.nanmean(auc[x]['early'][150:200])) for x in mouselist]
plt.bar(0.5,np.mean(data),width=1,color='grey')
plt.errorbar(0.5,np.mean(data), np.std(data)/np.sqrt(len(mouselist)),color='k',capsize=3,linewidth=0.5)
[plt.scatter(random.uniform(0,0.8)+0.1,x,1,'k') for x in data]
plt.plot([-0.5,1.5],[0,0],'k:',linewidth=0.35)
plt.xlim([-0.5,1.5])
plt.ylim([-0.6,0.6])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xticks([0.5],[])
plt.yticks([-0.6,-0.3,0,0.3,0.6])
plt.ylabel('\u0394Nrom. DA response\n(CS2-CS1)')
fig.savefig(dir+'//backpropagation_delta_2cues.pdf',bbox_inches='tight')

plt.close('all')
clr_light = ['grey','lightblue']
clr = ['k','b']
fig = plt.figure(figsize=(2.25*cm, 3*cm))
rect = 0.6,0.1,0.4,0.9
rect = [x*0.85 for x in rect]
ax = fig.add_axes(rect)
data = [np.divide(np.nanmean(auc[x]['late'][150:200]),np.nanmean(auc[x]['early'][150:200])) for x in mouselist]
plt.bar(0.5,np.mean(data),width=1,color='grey')
[plt.scatter(random.uniform(0,0.8)+0.1,x,1,'k') for x in data]
plt.errorbar(0.5,np.mean(data),np.std(data)/np.sqrt(len(mouselist)),color='k',capsize=3,linewidth=0.5)
plt.plot([-0.5,1.5],[0,0],'k:',linewidth=0.35)
plt.xlim([-0.5,1.5])
plt.ylim([-0.1,1])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xticks([0.5],[])
plt.yticks([0,0.5,1])
plt.gcf().set_size_inches(2.25*cm, 3*cm)
plt.ylabel('\u0394Norm. CS2 response\n(1 = CS1 response)')
fig.savefig(dir+'//backpropagation_2cues_last50_2.pdf',bbox_inches='tight')


