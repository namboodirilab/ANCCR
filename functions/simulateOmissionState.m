function [eventlog,cueoutcomelog] = simulateOmissionState(eventlog, cue_label, omission_label, rew_prob)
%SIMULATEOMISSIONSTATE: Output eventlog and cueoutcomelog for the given cue/reward  
%parameters cueoutcomelog is cell, which size is determined by
%number of cue identities. Each cell is consisted of eventlog index of cue,
%and outcome identity of every incidence of cue, whether it is followed by
%reward (1), recognized omission of reward (0), or unrecognized omission of
%reward (-1). Recognization of omission state is dependent on reward
%probability -- if reward probability is high, animals easily recognize
%omission of reward as well. 
cueoutcomelog = cell(length(cue_label),1);
for icue = 1:length(cue_label)
    if ~isnan(omission_label(icue))
        omidx = find(eventlog(:,1)==omission_label(icue));
        cueidx = find(eventlog(:,1)==cue_label(icue));
        cueoutcomelog{icue} = ones(length(cueidx),1);
        cueomission = find(cellfun(@(x) ismember(x+1,omidx),num2cell(cueidx)));
        cueoutcomelog{icue}(cueomission) = 0;        
        if ~isempty(omidx)
            outom = sort(randsample(length(omidx),round(length(omidx)*(1-rew_prob(icue)))));
            cueoutcomelog{icue}(cueomission(outom)) = -1;   
            eventlog(omidx(outom),:) = [];
        end
    end
end
for icue = 1:length(cue_label)
% 1st column: eventlog index
% 2nd column: outcome identity; 1, rewarded; 0, omission-recognized; -1, omission-not recognized
% For the calculation of omission response, need to have zero
% response in omission trials that are not recognized by an animal
cueoutcomelog{icue} = [find(eventlog(:,1)==cue_label(icue)),cueoutcomelog{icue}];
end
