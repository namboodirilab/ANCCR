function [eventlog,IRI] = simulateEvents(n_cues, cue_label, reward_label, ...
    reward_mag, omissionlabel, mean_ITI, max_ITI, cue_rew_delay, rew_prob,...
    postrewdelay, bgdrw_label, bgdrw_IRI, bgdrw_cue_delay, bgdrw_mag)
%SIMULATEEVENTS: Output an eventlog for the given cue reward parameters. 

% First check if one cue or multiple are being simulated (scalar vs. array)
if length(mean_ITI)==1
    mean_ITI = repmat(mean_ITI,1,length(cue_label));
end
if length(max_ITI)==1
    max_ITI = repmat(max_ITI,1,length(cue_label));
end
if length(cue_rew_delay)==1
    cue_rew_delay = repmat(cue_rew_delay,1,length(cue_label));
end
if length(postrewdelay)==1
    postrewdelay = repmat(postrewdelay,1,length(cue_label));
end
if length(reward_mag)==1
    reward_mag = repmat(reward_mag,1,length(cue_label));
end

eventlog = NaN(2*sum(n_cues),3);
order_cue = cell2mat(cellfun(@(x,y) ones(x,1)*y,...
    num2cell(n_cues(:)),num2cell([1:length(cue_label)]'),'UniformOutput',false)); % order of cues
order_cue = order_cue(randperm(length(order_cue)));

running_time = 0;
running_idx = 0;
% Loop case for keeping track of omitted trials
for i = 1:sum(n_cues)
    running_idx = running_idx + 1;
%     if length(unique(cue_label)) == length(cue_label)
%         icue = cue_label(order_cue(i));
%     else
        icue = order_cue(i);
%     end
    new_ts = exprnd(mean_ITI(icue));
    if new_ts > max_ITI(icue)
        new_ts = max_ITI(icue);
    end
    eventlog(running_idx, 1) = cue_label(icue);
    eventlog(running_idx, 2) = new_ts + running_time;
    eventlog(running_idx, 3) = 0;
    running_time = running_time + new_ts;
    if rew_prob(icue) > rand(1)
        running_idx = running_idx + 1;
        eventlog(running_idx, 1) = reward_label(icue);
        eventlog(running_idx, 2) = running_time + cue_rew_delay(icue);
        eventlog(running_idx, 3) = reward_mag(icue);
        running_time = running_time + cue_rew_delay(icue);
    else
        if ~isnan(omissionlabel(icue))
            running_idx = running_idx + 1;
            eventlog(running_idx, 1) = omissionlabel(icue);
            eventlog(running_idx, 2) = running_time + cue_rew_delay(icue);
            eventlog(running_idx, 3) = 0;
            running_time = running_time + cue_rew_delay(icue);
        end
    end
    running_time = running_time + postrewdelay(icue);
end

meanrewardnum = nan(1,length(cue_label));
for icue = 1:length(cue_label)
    meanrewardnum(icue) = rew_prob(icue); 
end

if nargin>10
    if ~isnan(bgdrw_IRI)
    numbgdrw = round(running_time/bgdrw_IRI)+100; % generate enough number of bgd reward 
    eventlog_bgdrw = simulateBackgroundRewards(numbgdrw, bgdrw_IRI, bgdrw_label, bgdrw_mag, 0);
    eventlog_bgdrw(eventlog_bgdrw(:,2)>running_time,:) = [];
    for icue = 1:length(cue_label)
        cuetimes = eventlog(eventlog(:,1)==cue_label(icue),2);
        in_bgdrw_cue_delay = sum(cell2mat(cellfun(@(x) x-eventlog_bgdrw(:,2)<=bgdrw_cue_delay &...
            x-eventlog_bgdrw(:,2)>=0,num2cell(cuetimes),'UniformOutput',false)'),2)>0;
        in_cue_rew_delay = sum(cell2mat(cellfun(@(x) eventlog_bgdrw(:,2)-x<=cue_rew_delay(icue)+postrewdelay(icue) &...
            eventlog_bgdrw(:,2)-x>=0,num2cell(cuetimes),'UniformOutput',false)'),2)>0;
        out = in_bgdrw_cue_delay | in_cue_rew_delay;
        eventlog_bgdrw(out,:) = [];
        meanrewardnum(icue) = meanrewardnum(icue)+(mean_ITI(icue)-bgdrw_cue_delay)/bgdrw_IRI;
    end
    eventlog = [eventlog;eventlog_bgdrw];
    eventlog = sortrows(eventlog,2);
    end
end

IRI = sum((cue_rew_delay+mean_ITI+postrewdelay).*(n_cues/n_cues(1)))/...
    sum(meanrewardnum.*(n_cues/n_cues(1)));

eventlog = rmmissing(eventlog);
end