function [time,alignedsignal] =...
    alignsignal2event(signaltime,signal,eventtime,window,resolution)
%ALIGNSIGNAL2EVENTS: align signal to each eventtime and generate aligned
%signal. Window constrained the start and end of aligned signal from each
%event time. If resolution is larger than zero, smooth alignedsignal using
%gaussian kernel with resolution. 

dt = round(nanmean(diff(signaltime))*10)/10; % find average bin size
winbin = round(window/dt); % window in unit of number of bin
time = [winbin(1):winbin(2)]*dt; 

% find the closest time bin from each event
[~,eventidxtemp] = cellfun(@(x) min(abs(signaltime-x)),num2cell(eventtime));

alignedsignal = NaN(length(eventtime),length(time));
out = isnan(eventtime); 
alignedsignaltemp = NaN(sum(~out),length(time));
eventidxtemp = eventidxtemp(~out);

% align signal
alignedsignaltemp(1,[winbin(1):winbin(2)]+eventidxtemp(1)>0) =...
    signal(max([1,eventidxtemp(1)+winbin(1)]):eventidxtemp(1)+winbin(2))';
outend = find(winbin(2)+eventidxtemp<length(signal),1,'last');
alignedsignaltemp(2:outend,:) = cell2mat(cellfun(@(x) signal(x+winbin(1):x+winbin(2))',...
    num2cell(eventidxtemp(2:outend)),'UniformOutput',false));
if outend<sum(~out)
    for iO = outend+1:sum(~out)
        alignedsignaltemp(iO,[winbin(1):winbin(2)]+eventidxtemp(iO)<=length(signal)) =...
            signal(eventidxtemp(iO)+winbin(1):min([length(signal),eventidxtemp(iO)+winbin(2)]))';
    end
end
alignedsignal(~out,:) = alignedsignaltemp;

% gaussian smoothing
if resolution>0
    alignedsignal = conv2(alignedsignal,fspecial('Gaussian',[1 5*resolution],resolution),'same');
end
end