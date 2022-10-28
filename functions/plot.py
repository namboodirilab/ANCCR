def align_events_to_reference(eventlog,eventindex,referenceindex,window,binsize,resolution):
    # align instances of one event to the other event
    # both of them are events that are stored in matlab eventlog

    # eventlog: eventlog from matlab data file
    # eventindex: the index of event you want to align
    # referenceindex
    #   - if referenceindex is an integer: the index of reference event you want to align to
    #   - if reference is an array: the time series of reference event you want to align to
    # window: the time window around each instance of reference event
    # binsize: bin size for histogram
    # resolution: determines the standard deviation for Gaussian kernel
    #   - sigma = resolution*binsize
    #   - if resolution is zero, not doing Gaussian convolution

    import numpy as np
    import scipy.ndimage as nd

    if isinstance(referenceindex, int):
        reference_times = eventlog[eventlog[:,0]==referenceindex,1]
    else:
        reference_times = referenceindex
    event_times = eventlog[eventlog[:,0]==eventindex,1]

    nbins = int(np.diff(window) / binsize)
    nrefs = len(reference_times)

    aligned_event = []
    hist_event = np.zeros(shape=(nrefs, nbins))
    for i, v in enumerate(reference_times):
        data = event_times[(event_times > v + window[0]) & (event_times < v + window[1])].tolist() - v
        aligned_event.append(data)
        hist_event[i, :] = np.histogram(data, nbins, range=(window[0], window[1]))[0] * (1000 / binsize)

    psthtime = np.transpose(range(window[0],window[1],binsize))

    meanpsth = np.mean(hist_event, 0)
    sempsth = np.std(hist_event, 0) / np.sqrt(hist_event.shape[0])
    if resolution > 0:
        meanpsth = nd.gaussian_filter1d(meanpsth, resolution)
        sempsth = nd.gaussian_filter1d(sempsth, resolution)

    return aligned_event, meanpsth, sempsth, psthtime

def align_signal_to_reference(signal,timestamp,eventlog,referenceindex,window,resolution,baselinenorm):
    # align continuous signal (like photometry trace) to the other event
    # this code assumes the signal was collected in a constant rate

    # signal: signal you want to align
    # timestamps: timestamps of signal
    # eventlog: eventlog from matlab data file
    # referenceindex
    #   - if referenceindex is an integer: the index of reference event you want to align to
    #   - if reference is an array: the time series of reference event you want to align to
    # window: the time window around each instance of reference event
    # resolution: determines the standard deviation for Gaussian kernel
    #   - sigma = resolution*(averaged interval b/w timestamps)
    #   - if resolution is zero, not doing Gaussian convolution
    # baselinenorm: whether normalize signals by the baseline, which is hard coded as [-1,0]

    import numpy as np
    import scipy.ndimage as nd

    if isinstance(referenceindex, int):
        reference_times = eventlog[eventlog[:, 0] == referenceindex, 1]
    else:
        reference_times = referenceindex

    baselinewindow = [-1000, 0]

    frameinterval = np.mean(np.diff(timestamp))
    framewindow = [int(i) for i in np.round(window/frameinterval)]
    framewindow_bl = [int(i) for i in np.round(baselinewindow / frameinterval)]
    nbins = np.sum(np.abs(framewindow))+1
    nrefs = len(reference_times)



    aligned_event = np.full((nrefs,nbins),np.nan)
    for i, v in enumerate(reference_times):
        closestframe = np.argmin(np.abs(timestamp-v))
        if framewindow[0]+closestframe<0:
            data = signal[np.arange(0,framewindow[1]+1+closestframe)]
            aligned_event[i,nbins-len(data):] = data
        elif framewindow[1]+1+closestframe>len(signal):
            data = signal[np.arange(framewindow[0]+closestframe,len(signal))]
            aligned_event[i,:-(nbins-len(data))] = data
        else:
            data = signal[np.arange(framewindow[0],framewindow[1]+1)+closestframe]
            aligned_event[i, :] = data
        if baselinenorm == 1:
            basemean = np.mean(signal[np.arange(framewindow_bl[0],framewindow_bl[1]+1)+closestframe])
            aligned_event[i, :] = data-basemean

    meanpsth = np.mean(aligned_event,0)
    sempsth = np.std(aligned_event,0)/np.sqrt(nrefs)
    if resolution > 0:
        meanpsth = nd.gaussian_filter1d(meanpsth, resolution)
        sempsth = nd.gaussian_filter1d(sempsth, resolution)
    psthtime = np.arange(framewindow[0],framewindow[1]+1)*frameinterval

    return aligned_event, meanpsth, sempsth, psthtime

def first_event_time_after_reference(eventlog,eventindex,referenceindex,window):
    # find time when event happens for the first time after each incidence of reference event
    # both of them are events that are stored in matlab eventlog

    # eventlog: eventlog from matlab data file
    # eventindex: the index of event you want to find
    # referenceindex
    #   - if referenceindex is an integer: the index of reference event
    #   - if reference is an array: the time series of reference event
    # window: searching window; if there is no event within this window around a reference event, have nan for output
    import numpy as np

    if isinstance(referenceindex, int):
        reference_times = eventlog[eventlog[:, 0] == referenceindex, 1]
    else:
        reference_times = referenceindex
    event_times = eventlog[eventlog[:,0]==eventindex,1]
    nrefs = len(reference_times)

    first_eventtime = np.full(nrefs,np.nan)
    for i,v in enumerate(reference_times):
        temp = [t for t in event_times if (t>v) and (t-v)<=window]
        if len(temp)>0:
            first_eventtime[i] = temp[0]
    return first_eventtime

def calculate_auc(signal,signaltime,reference_times,window):
    # calculate area under curve of signal during specified window from reference event

    # referenceindex
    #   - if referenceindex is an integer: the index of reference event
    #   - if reference is an array: the time series of reference event
    # window: window for calculation of auc
    import numpy as np

    binsize = np.mean(np.diff(signaltime))
    auc = [np.sum(signal[np.logical_and(signaltime >= v + window[0], signaltime <= v + window[1])])*binsize for v in reference_times]
    #auc = [np.trapz(signal[np.logical_and(signaltime >= v + window[0], signaltime <= v + window[1])],
    #                signaltime[np.logical_and(signaltime >= v + window[0], signaltime <= v + window[1])]) for v in reference_times]
    n = [np.sum(np.logical_and(signaltime >= v + window[0], signaltime <= v + window[1])) for v in reference_times]
    auc = [np.nan if x==0 else y for x,y in zip(n,auc)]
    return auc

def calculate_numevents(eventlog,eventindex,referenceindex,window):
    # calculate number of events during specified window from reference_event

    import numpy as np
    if isinstance(referenceindex, int):
        reference_times = eventlog[eventlog[:, 0] == referenceindex, 1]
    else:
        reference_times = referenceindex
    event_times = eventlog[eventlog[:,0]==eventindex,1]
    nevent = [np.sum(np.logical_and(event_times>=i+window[0],event_times<i+window[1])) for i in reference_times]
    return nevent


def plot_events(eventlog,eventindex,referenceindex,window,binsize,resolution,clr,ylabels,fig):
    import numpy as np
    (ax1, ax2) = fig.subplots(2, 1)
    itrial = 1
    for i,v in enumerate(referenceindex):
        aligned_event, meanpsth, sempsth, psthtime = align_events_to_reference(eventlog, eventindex, v, window, binsize, resolution)
        ax1.scatter(np.concatenate(aligned_event)/1000,np.concatenate([np.ones(np.shape(v),dtype=int)*(i+itrial) for i,v in enumerate(aligned_event)]),
                    c=clr[i],s=0.5,edgecolor=None)
        ax2.fill_between(psthtime/1000,meanpsth+sempsth,meanpsth-sempsth,alpha=0.3,facecolor=clr[i],linewidth=0)
        ax2.plot(psthtime/1000,meanpsth,color=clr[i])
        itrial = itrial+len(aligned_event)
    ax1.plot([0,0],[0,itrial],'k:',linewidth=0.35)
    ax2ylim = [np.floor(ax2.get_ylim()[0]),np.ceil(ax2.get_ylim()[1])]
    ax2.plot([0,0],ax2ylim,'k:',linewidth=0.35)
    ax1.set_xlim(window[0]/1000,window[1]/1000)
    ax2.set_xlim(window[0]/1000,window[1]/1000)
    ax1.set_ylim(0.5,itrial+0.5)
    ax2.set_ylim(ax2ylim[0],ax2ylim[1])
    ax2.set_xlabel('Time (s)')
    ax1.set_ylabel(ylabels[0])
    ax2.set_ylabel(ylabels[1])
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    return ax1, ax2


def plot_signal(signal,timestamp,eventlog,referenceindex,window,resolution,clr,ylabels,fig,baselinenorm):
    import numpy as np
    import matplotlib.pyplot as plt
    (ax1, ax2) = fig.subplots(2, 1)
    itrial = 1
    for i, v in enumerate(referenceindex):
        aligned_event, meanpsth, sempsth, psthtime = align_signal_to_reference(signal,timestamp,eventlog,v,window,resolution,baselinenorm)
        dx = (psthtime[1] - psthtime[0]) / 2
        im = ax1.imshow(aligned_event, extent = [(window[0]-dx)/1000, (window[-1]+dx)/1000, itrial-0.5, np.shape(aligned_event)[0]+itrial+0.5],
                   aspect='auto')
        ax2.fill_between(psthtime / 1000, meanpsth + sempsth, meanpsth - sempsth, alpha=0.3, facecolor=clr[i], linewidth=0)
        ax2.plot(psthtime / 1000, meanpsth, color=clr[i])
        itrial = itrial + len(aligned_event)
    ax1.plot([0, 0], [0, itrial], 'k:',linewidth=0.35)
    ax2ylim = [np.floor(ax2.get_ylim()[0]), np.ceil(ax2.get_ylim()[1])]
    ax2.plot([0, 0], ax2ylim, 'k:',linewidth=0.35)
    ax1.set_xlim(window[0] / 1000, window[1] / 1000)
    ax2.set_xlim(window[0] / 1000, window[1] / 1000)
    ax1.set_ylim(0.5, itrial + 0.5)
    ax2.set_ylim(ax2ylim[0], ax2ylim[1])
    ax2.set_xlabel('Time (s)')
    ax1.set_ylabel(ylabels[0])
    ax2.set_ylabel(ylabels[1])
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    return ax1, ax2, im


def cumulative_analysis(signal,xmaxwrtchangetrial):
    import numpy as np
    from numpy.linalg import norm
    signalx = np.arange(1,len(signal)+1)/(len(signal)+1)
    cumsignal = np.cumsum(signal)/np.sum(signal)
    d = [norm(np.cross(np.subtract([1,1],[0,0]),np.subtract([x,y],[1,1])))/norm(np.subtract([1,1],[0,0])) for x,y in zip(signalx,cumsignal)]
    abruptness = max(d)
    changetrial = d.index(abruptness)

    if ~np.isnan(xmaxwrtchangetrial):
        cumsignal = np.cumsum(signal[:round(changetrial*xmaxwrtchangetrial)])/np.sum(signal[:round(changetrial*xmaxwrtchangetrial)])
        signalx = np.arange(1,round(changetrial*xmaxwrtchangetrial)+1)/(round(changetrial*xmaxwrtchangetrial)+1)
        d = [norm(np.cross(np.subtract([1,1],[0,0]),np.subtract([x,y],[1,1])))/norm(np.subtract([1,1],[0,0])) for x,y in zip(signalx,cumsignal)]
        abruptness = max(d)

    return abruptness, changetrial, cumsignal, signalx
