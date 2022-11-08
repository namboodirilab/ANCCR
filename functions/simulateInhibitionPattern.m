function [optolog] = simulateInhibitionPattern(eventlog, ...
    inhib_targets, inhib_mcts, inhib_mag, trials_off)
%SIMULATEINHIBITION generates 'optolog', an array where the first column
%includes indicators for whether inhibition is occuring for a given
%timestep and the second column encodes the magnitude of the 'substituted'
%DA response.

optolog = zeros(length(eventlog), 2);
assert(length(inhib_targets) == length(inhib_mag))
% If you only want to inhibit targets that follow a particular event/mct
if ~isnan(inhib_mcts)
    for i = 1:length(inhib_targets)
        inhib_targ_idxs = squeeze(find(eventlog(:,1) == inhib_targets(i))); % List of all rewards
        inhib_mct_idxs = squeeze(find(eventlog(:,1) == inhib_mcts(i))); % List of all C1s
        targ_idx = nan(length(inhib_targ_idxs),1); % This is meant to be a list of the indices that should be inhibited
        for j = 1:length(inhib_targ_idxs) % Loop through rewards and save the ones that are 'caused' by C1
            cause_idx = find((inhib_targ_idxs(j) - inhib_mct_idxs) == 1, 1);
            if isempty(cause_idx) % If a cause can't be found, don't save target
                targ_idx(j) = nan;
            else
                % Save index if cause is found
                targ_idx(j) = inhib_targ_idxs(j);
            end
        end
        targ_idx = rmmissing(targ_idx); % Remove nan for 'not found' instances
        count_to_trials_off = 0;
        for k = 1:length(targ_idx)
            count_to_trials_off = count_to_trials_off + 1;
            % Once intended trial num is reached, inhibit
            if count_to_trials_off == trials_off
                count_to_trials_off = 0;
                optolog(targ_idx(k),2) = inhib_mag(i);
                optolog(targ_idx(k),1) = 1;
                optolog(targ_idx(k)-1,1) = 2; 
            end
        end
    end
% If you want to generally inhibit all instances of an event
else
    for i = 1:length(inhib_targets)
        inhib_targ_idxs = squeeze(find(eventlog(:,1) == inhib_targets(i)));
        for j = 1:length(inhib_targ_idxs)
            optolog(inhib_targ_idxs(j), 2) = inhib_mag(i);
            optolog(inhib_targ_idxs(j), 1) = 1;
        end
    end
end
end