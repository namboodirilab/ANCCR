function [eventlog] = simulateEventChain(n_sequences, n_cues, omissionlabel, mean_ITI, ...
    max_ITI, intercuedelay, cuerewdelay, rew_prob, postrewdelay)
%SIMULATEEVENTCHAIN: Output an eventlog for the given cue/reward params.
eventlog = NaN(n_sequences*n_cues,2);
running_time = 0;
running_idx = 0;
for i = 1:n_sequences
    running_idx = running_idx + 1;
    new_ts = exprnd(mean_ITI);
    if new_ts > max_ITI
        new_ts = max_ITI;
    end
    running_time = running_time + new_ts;
    eventlog(running_idx, 1) = 1;
    eventlog(running_idx, 2) = running_time;
    eventlog(running_idx, 3) = 0;
    % Add all following cues
    for j = 2:n_cues
        running_idx = running_idx + 1;
        running_time = running_time + intercuedelay;
        eventlog(running_idx, 1) = j;
        eventlog(running_idx, 2) = running_time;
        eventlog(running_idx, 3) = 0;
    end
    % And add reward
    running_time = running_time + cuerewdelay;
    if rew_prob > rand(1)
        running_idx = running_idx + 1;
        eventlog(running_idx, 1) = j + 1;
        eventlog(running_idx, 2) = running_time;
        eventlog(running_idx, 3) = 1;
    else
        if ~isnan(omissionlabel)
            running_idx = running_idx + 1;
            eventlog(running_idx, 1) = omissionlabel;
            eventlog(running_idx, 2) = running_time;
            eventlog(running_idx, 3) = 0;
        end
    end
    running_time = running_time + postrewdelay;
end
eventlog = rmmissing(eventlog);
end