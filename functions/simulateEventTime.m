function [eventlog] = simulateEventTime(sim_time, cue_label, reward_label, ...
    rew_mag, mean_ITI, max_ITI, cue_rew_delay, rew_prob, post_rew_delay)
%SIMULATEEVENTTIME: Output an eventlog and omission indices for the given cue/
%reward parameters. Simulation runs for specified time rather than
%specific number of events.

eventlog = NaN(5000,3); % Picked arbitrarily large number for size, modify if too small
running_time = 0;
running_idx = 0;
while running_time < sim_time
    running_idx = running_idx + 1;
    new_ts = exprnd(mean_ITI);
    if new_ts > max_ITI
        new_ts = max_ITI;
    end
    running_time = running_time + new_ts; % cue delivery time
    eventlog(running_idx, 1) = cue_label;
    eventlog(running_idx, 2) = running_time;
    eventlog(running_idx, 3) = 0;

    running_time = running_time + cue_rew_delay; % reward delivery / omission time
    if rew_prob > rand(1)
        running_idx = running_idx + 1;
        eventlog(running_idx, 1) = reward_label;
        eventlog(running_idx, 2) = running_time;
        eventlog(running_idx, 3) = rew_mag;
    end

    running_time = running_time + post_rew_delay; % trial end time - next trial starts from here
end
eventlog = rmmissing(eventlog);
end