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

eventlog = NaN(2*sum(n_cues),3); % Initialize eventlog
% For multiple unique cues, generate random sequence to present cues
order_cue = cell2mat(cellfun(@(x,y) ones(x,1)*y,...
    num2cell(n_cues(:)),num2cell([1:length(cue_label)]'),'UniformOutput',false));
order_cue = order_cue(randperm(length(order_cue)));

running_time = 0;
running_idx = 0;
for i = 1:sum(n_cues)
    running_idx = running_idx + 1;
    icue = order_cue(i); % Index of unique cue being generated
    new_ts = exprnd(mean_ITI(icue)); % Exp. dist. timestamp for given cue
    if new_ts > max_ITI(icue)
        % If new ts exceeds max, set to max
        new_ts = max_ITI(icue);
    end
    % Update eventlog with cue label and new timestamp
    eventlog(running_idx, 1) = cue_label(icue);
    eventlog(running_idx, 2) = new_ts + running_time;
    eventlog(running_idx, 3) = 0;
    running_time = running_time + new_ts; % Update running time
    if rew_prob(icue) > rand(1) % Rng to determine if reward is delivered
        running_idx = running_idx + 1;
        eventlog(running_idx, 1) = reward_label(icue);
        eventlog(running_idx, 2) = running_time + cue_rew_delay(icue);
        eventlog(running_idx, 3) = reward_mag(icue);
        running_time = running_time + cue_rew_delay(icue);
    else
        % If no reward is delivered and omission is being tracked, save 
        % omission instance in eventlog
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

% Estimate number of rewards based on reward probability
meanrewardnum = nan(1,length(cue_label));
for icue = 1:length(cue_label)
    meanrewardnum(icue) = rew_prob(icue); 
end

% If background rewards are also being simulated...
if nargin>10
    if ~isnan(bgdrw_IRI)
        % Generate sufficient backgoround rewards (based on running time)
        numbgdrw = round(running_time/bgdrw_IRI)+100;
        eventlog_bgdrw = simulateBackgroundRewards(numbgdrw, bgdrw_IRI, bgdrw_label, bgdrw_mag, 0);
        eventlog_bgdrw(eventlog_bgdrw(:,2)>running_time,:) = [];
        for icue = 1:length(cue_label)
            cuetimes = eventlog(eventlog(:,1)==cue_label(icue),2);
            % Find instances where background rewards are too close to cues
            in_bgdrw_cue_delay = sum(cell2mat(cellfun(@(x) x-eventlog_bgdrw(:,2)<=bgdrw_cue_delay &...
                x-eventlog_bgdrw(:,2)>=0,num2cell(cuetimes),'UniformOutput',false)'),2)>0;
            % Find instances where background rewards are in cue/rew delay
            in_cue_rew_delay = sum(cell2mat(cellfun(@(x) eventlog_bgdrw(:,2)-x<=cue_rew_delay(icue)+postrewdelay(icue) &...
                eventlog_bgdrw(:,2)-x>=0,num2cell(cuetimes),'UniformOutput',false)'),2)>0;
            out = in_bgdrw_cue_delay | in_cue_rew_delay;
            % Remove instances of background rewards that satisfy above
            eventlog_bgdrw(out,:) = [];
            meanrewardnum(icue) = meanrewardnum(icue)+(mean_ITI(icue)-bgdrw_cue_delay)/bgdrw_IRI;
        end
        % Update eventlog with background rewards and sort by time
        eventlog = [eventlog;eventlog_bgdrw];
        eventlog = sortrows(eventlog,2);
    end
end

% Calculate the mean interreward interval
IRI = sum((cue_rew_delay+mean_ITI+postrewdelay).*(n_cues/n_cues(1)))/...
    sum(meanrewardnum.*(n_cues/n_cues(1)));
% Remove nans from final eventlog
eventlog = rmmissing(eventlog);
end