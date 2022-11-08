function auc = aucsignal2event(signaltime,signal,eventtime,window,out)
%AUCSIGNAL2EVENT: calculate area under curve (AUC) of signal during window 
% around each event time. 

if nargin<5
    out = false(length(eventtime),1);
end

dt = nanmean(diff(signaltime)); % find average frame size
winbin = round(window/dt); % window in unit of number of bin

out = isnan(eventtime) | out;

% find the closest time bin from each event
[~,eventidxtmp] = cellfun(@(x) min(abs(signaltime-x)),num2cell(eventtime));
eventidxtmp(out) = NaN;
out(eventidxtmp+winbin(2)>length(signaltime)) = true;
out(eventidxtmp+winbin(1)<0) = true;

% calculate AUC
auc = NaN(length(eventtime),1);
auc(~out) = cellfun(@(x) sum(signal(x+winbin(1):x+winbin(2))),num2cell(eventidxtmp(~out)));

end