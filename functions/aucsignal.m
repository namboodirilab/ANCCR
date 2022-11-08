function auc = aucsignal(signal,time,window,out)
%AUCSIGNAL: calculate area under curve (AUC) of signal during window.
%Signal is a matrix with a dimension of observation X time, so this
%function iterate calculation of AUC for every observations.

if nargin<4
    out = false(size(signal,1),1);
end
binsize = nanmean(diff(time)); % find average frame size
auc = NaN(size(signal,1),1);
auc(~out) = nansum(signal(~out,time>=window(1) & time<=window(2)),2)*binsize;
end

