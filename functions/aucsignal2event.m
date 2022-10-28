function auc = aucsignal2event(signaltime,signal,eventtime,window,out)
if nargin<5
    out = false(length(eventtime),1);
end

framesize = nanmean(diff(signaltime));
winframe = round(window/framesize);

out = isnan(eventtime) | out;
[~,eventidxtmp] = cellfun(@(x) min(abs(signaltime-x)),num2cell(eventtime));
eventidxtmp(out) = NaN;
out(eventidxtmp+winframe(2)>length(signaltime)) = true;
out(eventidxtmp+winframe(1)<0) = true;

auc = NaN(length(eventtime),1);
auc(~out) = cellfun(@(x) sum(signal(x+winframe(1):x+winframe(2))),num2cell(eventidxtmp(~out)));

end