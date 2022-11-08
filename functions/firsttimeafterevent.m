function firsttime = firsttimeafterevent(targeteventtime,referenceeventtime,window)
% FIRSTTIMEAFTEREVENT: find first target event time after each reference
% event. Window specifies the searching window from referenceeventtime.

if isnan(window)
    % if window is nan, find time of first target event after reference
    % event time w/o limited window 
    firsttime = cellfun(@(x) targeteventtime(find(targeteventtime>=x,1,'first')),...
        num2cell(referenceeventtime),'UniformOutput',false);
else
    if length(window)==length(referenceeventtime)
        % if length of window is length of referenceeventtime, search target
        % event until window of each referenceeventtime
        firsttime = cellfun(@(x,y) targeteventtime(find(targeteventtime>x & targeteventtime<y,1,'first')),...
            num2cell(referenceeventtime),num2cell(window),'UniformOutput',false);
    elseif length(window)==1
        % if window is a single value, search target event during [0, window]
        % from each referenceeventtime
        firsttime = cellfun(@(x) targeteventtime(find(targeteventtime>x & targeteventtime<x+window,1,'first')),...
            num2cell(referenceeventtime),'UniformOutput',false);
    end
end
% if target event doesn't happen within window, have nan
firsttime(cellfun(@isempty,firsttime)) = {NaN}; 
firsttime = cell2mat(firsttime);
