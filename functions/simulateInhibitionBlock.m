function [optolog, first_inh] = simulateInhibitionBlock(eventlog, inhib_targets, inhib_mcts, inhib_mag, inhib_prob, start, stop)
%SIMULATEINHIBITION generates 'optolog', an array where the first column
%includes indicators for whether inhibition is occuring for a given
%timestep and the second column encodes the magnitude of the 'substituted'
%DA response. 
optolog = zeros(length(eventlog), 2);
assert(length(inhib_targets) == length(inhib_mag))

if ~isnan(inhib_mcts) % If inhibiting rewards after MCTs
    for i = 1:length(inhib_targets)
        inhib_targ_idxs = find(eventlog(:,1) == inhib_targets(i));
        if isnan(stop)
            stop = length(inhib_targ_idxs);
        end
        inhib_mct_idxs = find(eventlog(:,1) == inhib_mcts(i));
        % Only inhibit after start idx
        inhib_targ_idxs = squeeze(inhib_targ_idxs(start:stop)); 
        inhib_mct_idxs = squeeze(inhib_mct_idxs(start:stop)); 
        for j = 1:length(inhib_targ_idxs)
            cause_idx = find((inhib_mct_idxs - inhib_targ_idxs(j)) == 1, 1);
            if ~isempty(cause_idx)
                dice = rand(); % Rng to determine if inhibition happens
                if dice <= inhib_prob
                    optolog(inhib_targ_idxs(j), 2) = inhib_mag(i);
                    % Inhibited targets get label 1
                    optolog(inhib_targ_idxs(j), 1) = 1;
                    % If a 'caused' target is successfully inhibited,
                    % the cause gets label 2
                    optolog(inhib_mct_idxs(cause_idx), 1) = 2;
                else
                    optolog(inhib_targ_idxs(j), 2) = 0;
                    % Uninhibited targets get label 1
                    optolog(inhib_targ_idxs(j), 1) = 0;
                    % If a target is caused but not successfully inhibited,
                    % the cause gets label 3
                    optolog(inhib_mct_idxs(cause_idx), 1) = 3;
                end
            else % If no causes are found, do not inhibit
                optolog(inhib_targ_idxs(j), 2) = 0;
                optolog(inhib_targ_idxs(j), 1) = 0;
            end
        end
    end
else % If generally inhibiting rew.
    for i = 1:length(inhib_targets)
        inhib_targ_idxs = find(eventlog(:,1) == inhib_targets(i));
        if isnan(stop)
            stop = length(inhib_targ_idxs);
        end
        inhib_targ_idxs = squeeze(inhib_targ_idxs(start:stop)); 
        for j = 1:length(inhib_targ_idxs)
            dice = rand(); % Rng to determine if inhibition happens
            if dice <= inhib_prob
                optolog(inhib_targ_idxs(j), 2) = inhib_mag(i);
                optolog(inhib_targ_idxs(j), 1) = 1;
            else
                optolog(inhib_targ_idxs(j), 2) = 0;
                optolog(inhib_targ_idxs(j), 1) = 0;
            end
        end
    end
end
    first_inh = find(optolog(:, 1), 1); % Save the first inhibition idx
end