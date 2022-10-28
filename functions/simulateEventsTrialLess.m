function [eventlog] = simulateEventsTrialLess(n_cues, cue_label, reward_label, ...
    reward_mag, omissionlabel, mean_ITI, max_ITI, min_ITI, cue_rew_delay, rew_prob)
%SIMULATEEVENTS: Output an eventlog and omission indices for the given cue/
%reward parameters.

if length(reward_mag)==1
    reward_mag = repmat(reward_mag,1,length(cue_label));
end
if length(omissionlabel)==1
    omissionlabel = repmat(omissionlabel,1,length(cue_label));
end
if length(mean_ITI)==1
    mean_ITI = repmat(mean_ITI,1,length(cue_label));
end
if length(max_ITI)==1
    max_ITI = repmat(max_ITI,1,length(cue_label));
end
if length(min_ITI)==1
    min_ITI = repmat(min_ITI,1,length(cue_label));
end

eventlog = NaN(2*sum(n_cues),3);

running_idx = 0;
for jcue = 1:length(cue_label)
    running_time = 0;
    for i = 1:n_cues(jcue)
        running_idx = running_idx + 1;
        new_ts = max_ITI(jcue)+1;
        while new_ts > max_ITI(jcue) | new_ts<min_ITI(jcue)
            new_ts = exprnd(mean_ITI(jcue));
        end
        eventlog(running_idx, 1) = cue_label(jcue);
        eventlog(running_idx, 2) = new_ts + running_time;
        eventlog(running_idx, 3) = 0;
        
        if rew_prob(jcue) > rand(1)
            running_idx = running_idx + 1;
            eventlog(running_idx, 1) = reward_label(jcue);
            eventlog(running_idx, 2) = new_ts + running_time + cue_rew_delay(jcue);
            eventlog(running_idx, 3) = reward_mag(jcue);
        else
            if ~isnan(omissionlabel(jcue))
                running_idx = running_idx + 1;
                eventlog(running_idx, 1) = omissionlabel(jcue);
                eventlog(running_idx, 2) = new_ts + running_time + cue_rew_delay(jcue);
                eventlog(running_idx, 3) = 0;
            end
        end
        running_time = running_time + new_ts;
    end
end
eventlog = rmmissing(eventlog);
eventlog = sortrows(eventlog,2);

% finish the session if the number of any cue reaches to n_cue
lasttrial = cellfun(@(x) find(eventlog(:,1)==x,1,'last'),num2cell(reward_label(rew_prob>0)));
eventlog(min(lasttrial)+1:end,:) = [];
end