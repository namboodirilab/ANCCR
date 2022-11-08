import functions.load as fnl
import functions.plot as fnp
import matplotlib.pyplot as plt
import numpy as np
import random

# name spaces for namlab nwb extension
namespaces = ["C:\\Users\\Huijeong Jeong\\ndx-photometry-namlab\\spec\\ndx-photometry-namlab.namespace.yaml",
              "C:\\Users\\Huijeong Jeong\\ndx-eventlog-namlab\\spec\\ndx-eventlog-namlab.namespace.yaml"]
# DANDI set id
dandiset_id = '000351'

# animal list
dathetlist = ['M2','M3','M4','M5','M6','M7']
datwtlist = ['F1']
wtlist = ['F1','F2','F3','M1','M2','M3']
mouselist = ['HJ-FP-datWT-stGtACR-'+i for i in datwtlist] + ['HJ-FP-WT-stGtACR-'+i for i in wtlist]

# event indices
csindex = 15
rewardindex = 10
lickindex = 5

auc = {}
for im, mousename in enumerate(mouselist):
    print(mousename)

    # load DANDI url of an animal
    url, path = fnl.load_dandi_url(dandiset_id, mousename)
    url= [y for x,y in zip(path,url) if 'Pavlovian' in x]

    auc[mousename] = {}
    auc[mousename]['early'] = [] # dopamine response to CS1
    auc[mousename]['late'] = [] # dopamine response to CS2
    for i, v in enumerate(url):
        # load eventlog and dff from nwb file
        results,_ = fnl.load_nwb(v, namespaces,[('a','eventlog'),('p','photometry','dff')])

        # CS timestamps
        cstime = results['eventlog']['eventtime'][results['eventlog']['eventindex']==csindex]
        cs1time = cstime[range(0,len(cstime),2)]
        cs2time = cstime[range(1,len(cstime),2)]

        # calculate AUC during CS1(early) and CS2(late)
        early_temp = fnp.calculate_auc(results['dff']['data'], results['dff']['timestamps'], cs1time, [0,1])
        late_temp = fnp.calculate_auc(results['dff']['data'], results['dff']['timestamps'],cs2time,[0,1])

        # in first session, calculate dopamine response to reward (first lick after reward delivery)
        if i==0:
            firstlicktimes = fnp.first_event_time_after_reference(results['eventlog']['eventindex'],results['eventlog']['eventtime'], lickindex, rewardindex, 5)
            auc[mousename]['reward'] = fnp.calculate_auc(results['dff']['data'], results['dff']['timestamps'], firstlicktimes, [0, 1])
        auc[mousename]['early'] = auc[mousename]['early'] + early_temp
        auc[mousename]['late'] = auc[mousename]['late'] + late_temp

## set plotting parameters
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

dir = 'D:\OneDrive - University of California, San Francisco\\figures\manuscript\dopamine_contingency\\revision\\test'
cm = 1/2.54

# fig 6D left
ntrials = 200
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

# fig 6D middle
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

# fig 6D right
clr_light = ['grey','lightblue']
clr = ['k','b']
fig = plt.figure(figsize=(2.25*cm, 3*cm))
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


