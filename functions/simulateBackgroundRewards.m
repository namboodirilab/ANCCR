function [eventlog] = simulateBackgroundRewards(numrewards, rewITI, rewlabel, rewmag, truncation)
%SIMULATEBACKGROUNDREWARDS outputs an eventlog with exclusively background
%rewards. Inputs may be scalar or arrays depending on the number of unique
%background rewards specified.
    % Set upper limit on ITI (only used if truncation == 1)
    maxrewITI = 3*rewITI;
    % Initialize empty eventlog
    eventlog = NaN(sum(numrewards),3);
    % Keep track of reward indices
    running_idx = 0; 
    % If multiple unique rewards, loop instance for each type
    for irw = 1:length(rewlabel)
        % Start simulation time at 0
        running_time = 0; 
        % For each reward (irw) simulate 'numrewards' timesteps
        for i = 1:numrewards(irw) 
            running_idx = running_idx + 1;
            % Generate exp. dist. timestep for current reward instance
            new_ts = exprnd(rewITI(irw)); 
            if truncation==1 % If capping ITI at maxITI
                while new_ts > maxrewITI(irw)
                    % Simulate a new ts if maxITI is exceeded
                    new_ts = exprnd(rewITI(irw));
                end
            end
            % Update eventlog
            eventlog(running_idx, 1) = rewlabel(irw); % Column 1 = label
            eventlog(running_idx, 2) = new_ts + running_time; % " 2 = ts
            eventlog(running_idx, 3) = rewmag(irw); % " 3 = rew. magnitude
            running_time = running_time + new_ts; % Update running time
        end
    end
    
    % Resort eventlog by timestamp
    eventlog = sortrows(eventlog,2);
    
    % Once any reward reaches to numrewards, finish the session
    lasttrial = min(cellfun(@(x) find(eventlog(:,1)==x,1,'last'),num2cell(rewlabel)));
    eventlog(lasttrial+1:end,:) = [];
end