function [eventlog] = simulateEventChain(n_sequences, n_cues, omissionlabel, mean_ITI, ...
    max_ITI, intercuedelay, cuerewdelay, rew_prob, postrewdelay)
%SIMULATEEVENTCHAIN: Output an eventlog for the given cue/reward params.
%Expects a sequence of cues preceding a reward.

eventlog = NaN(n_sequences*n_cues,2); % Initialize eventlog
running_time = 0; % Initialize time
running_idx = 0; % Initialize index counter
for i = 1:n_sequences % Loop through unique sequences
    running_idx = running_idx + 1; % Increment running idx
    new_ts = exprnd(mean_ITI); % Generate ts from exp. distribution
    if new_ts > max_ITI
        % If ts exceeds maxITI, set ts to maxITI
        new_ts = max_ITI;
    end
    running_time = running_time + new_ts; % Update running time
    % Update simulation log
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
        % If keeping track of omission, save omission instances
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