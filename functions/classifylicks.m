function [lickboutidx,lickconsumidx,consumboutlength,lickrwidx] =...
    classifylicks(licktime,rewardtime,sessionendtime,boutinterval)

% classify licks into bout onset, bout offset, and mid-bout licks.
% consumption bout is defined as the lick bout containing the first lick
% after reward delivery. If consumption bout has started before reward
% delivery, count only licks after reward delivery 
% lickrwidx indicates the latest reward number from current lick

firstlicktime = firsttimeafterevent(licktime,rewardtime,[rewardtime(2:end);sessionendtime]);
try
rewardwin = [[0;rewardtime],[rewardtime(1:end);sessionendtime]];
nReward = length(rewardtime);
lickrwidx = sum(cell2mat(cellfun(@(x,y) double(licktime>=x(1) & licktime<x(2))*y,...
    mat2cell(rewardwin,ones(nReward+1,1),2),num2cell([1:nReward+1]'),'UniformOutput',false)'),2)-1;
catch
nReward
end
lickboutidx = NaN(length(licktime),1);
lickdiff = diff([0; licktime]);
lickboutidx(lickdiff>=boutinterval) = 1; %bout onset lick
lickboutidx([lickdiff(2:end);0]>=boutinterval) = 3; %bout offset lick
lickboutidx(isnan(lickboutidx)) = 2; %mid-bout lick

nLick = length(licktime);
lickconsumidx = false(nLick,1);
consumboutlength = NaN(nReward,1);
for iR = 1:nReward
    firstlickidx = find(licktime==firstlicktime(iR));
    if isempty(firstlickidx)
        continue;
    end
    lickboutend = find(licktime>firstlicktime(iR) & lickboutidx==3,1,'first');
    lickconsumidx(firstlickidx:lickboutend) = true;

    % for the last lick bout
    if sum(lickboutidx(licktime>firstlicktime(iR))==2)==sum(licktime>firstlicktime(iR))
        lickconsumidx(firstlickidx:end) = true;
    end
    
    % if consumption bout doesn't end until next reward delivery, calculate
    % the length of consumption bout of current reward as the last lick
    % before next reward and the first lick after current reward
    if iR==nReward
        consumboutendtime = min([sessionendtime,...
            licktime(find(lickboutidx==3&lickrwidx>=iR,1,'first'))]);
    else
        consumboutendtime = min([rewardtime(iR+1),sessionendtime,...
            licktime(find(lickboutidx==3&lickrwidx>=iR,1,'first'))]);
    end
    consumboutlength(iR) = consumboutendtime-firstlicktime(iR);
    lickconsumidx(licktime>=firstlicktime(iR) & licktime<=consumboutendtime) = 1;
end
end