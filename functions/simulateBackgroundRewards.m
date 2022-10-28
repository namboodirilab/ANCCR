function [eventlog] = simulateBackgroundRewards(numrewards, rewITI, rewlabel, rewmag, truncation)
%SIMULATEBACKGROUNDREWARDS
    maxrewITI = 3*rewITI;
    eventlog = NaN(sum(numrewards),3);
    
    running_idx = 0;
    for irw = 1:length(rewlabel)
        running_time = 0;
        for i = 1:numrewards(irw)
            running_idx = running_idx + 1;
            new_ts = exprnd(rewITI(irw));
            if truncation==1
                if new_ts > maxrewITI(irw)
                    new_ts = exprnd(rewITI(irw));
                end
            end
            eventlog(running_idx, 1) = rewlabel(irw);
            eventlog(running_idx, 2) = new_ts + running_time;
            eventlog(running_idx, 3) = rewmag(irw);
            running_time = running_time + new_ts;
        end
    end
    
    eventlog = sortrows(eventlog,2);
    
    % if one reward reaches to its numrewards, finish the session
    lasttrial = min(cellfun(@(x) find(eventlog(:,1)==x,1,'last'),num2cell(rewlabel)));
    eventlog(lasttrial+1:end,:) = [];
end