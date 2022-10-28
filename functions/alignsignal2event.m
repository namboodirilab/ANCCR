function [time,alignedsignal] =...
    alignsignal2event(signaltime,signal,eventtime,window,resolution)
framesize = round(nanmean(diff(signaltime))*10)/10;
winframe = round(window/framesize);
time = [winframe(1):winframe(2)]*framesize;

[~,eventidxtemp] = cellfun(@(x) min(abs(signaltime-x)),num2cell(eventtime));

alignedsignal = NaN(length(eventtime),length(time));
out = isnan(eventtime);

alignedsignaltemp = NaN(sum(~out),length(time));
eventidxtemp = eventidxtemp(~out);
alignedsignaltemp(1,[winframe(1):winframe(2)]+eventidxtemp(1)>0) =...
    signal(max([1,eventidxtemp(1)+winframe(1)]):eventidxtemp(1)+winframe(2))';
outend = find(winframe(2)+eventidxtemp<length(signal),1,'last');
alignedsignaltemp(2:outend,:) = cell2mat(cellfun(@(x) signal(x+winframe(1):x+winframe(2))',...
    num2cell(eventidxtemp(2:outend)),'UniformOutput',false));
if outend<sum(~out)
    for iO = outend+1:sum(~out)
        alignedsignaltemp(iO,[winframe(1):winframe(2)]+eventidxtemp(iO)<=length(signal)) =...
            signal(eventidxtemp(iO)+winframe(1):min([length(signal),eventidxtemp(iO)+winframe(2)]))';
    end
end
alignedsignal(~out,:) = alignedsignaltemp;

if resolution>0
alignedsignal = conv2(alignedsignal,fspecial('Gaussian',[1 5*resolution],resolution),'same');
end
end