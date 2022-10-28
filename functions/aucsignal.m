function auc = aucsignal(signal,time,window,out)
if nargin<4
    out = false(size(signal,1),1);
end
% out = out | sum(isnan(signal(:,time>=window(1) & time<=window(2))),2)>0;
binsize = nanmean(diff(time));
auc = NaN(size(signal,1),1);
auc(~out) = nansum(signal(~out,time>=window(1) & time<=window(2)),2)*binsize;
end

