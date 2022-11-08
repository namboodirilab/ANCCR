function [sessionendtime,licktime,CS,CStime,CSduration,fxrwtime,fxreward,...
    bgdrwtime,firstfxlicktime] = eventfrompavlovian(eventtime,eventindex,nosolenoidflag,trialless)
% calculate task related information from pavlovian task

if nargin<4
    trialless=false;
end

% event indices
index_lick = 5;
index_fxrw = 10;
index_bgdrw = 7;
index_CS1 = 15;
index_CS2 = 16;
index_sessionend = 0;
index_trialend = 14;

sessionendtime = eventtime(eventindex==index_sessionend); % session end time
if isempty(sessionendtime)
    sessionendtime = eventtime(end)+5000;
end

licktime = eventtime(eventindex==index_lick); % lick time

[CSidx,CStmp] = ismember(eventindex,[index_CS1 index_CS2]);
CS = CStmp(CSidx); % CS identity of each trial
CStime = eventtime(CSidx); % CS time

fxrwtime = eventtime(eventindex==index_fxrw); % cue-associated reward time
fxreward = 1-nosolenoidflag(eventindex==index_fxrw); % to account for reward omisson trials
bgdrwtime = eventtime(eventindex==index_bgdrw); % background reward time

CS(CStime>fxrwtime(end)) = [];
CStime(CStime>fxrwtime(end)) = [];
CSduration = nanmean(fxrwtime-CStime)-1000; % CS duration: 1s shorter than cue-reward delay; this is hard-coded 

trialendtime = eventtime(eventindex==index_trialend); % trial end time

% find first lick time after cue-associated reward delivery
if ~trialless
    if ~isempty(trialendtime)
        if nanmean(trialendtime-fxrwtime)>=2000
            % if trialendtime was logged & delay b/w trialendtime and reward
            % time was larger than 2s, search first lick until trialendtime
            firstfxlicktime = firsttimeafterevent(licktime,fxrwtime,trialendtime);
        else
            % else, search first lick until 1s before next cue delivery
            firstfxlicktime = firsttimeafterevent(licktime,fxrwtime,[CStime(2:end)-1000;sessionendtime]);
        end
    else
        firstfxlicktime = firsttimeafterevent(licktime,fxrwtime,[CStime(2:end)-1000;sessionendtime]);
    end
else
    % in trial-less task, find first lick after reward time w/o any
    % constrained window
    firstfxlicktime = firsttimeafterevent(licktime,fxrwtime,nan);
end
end
