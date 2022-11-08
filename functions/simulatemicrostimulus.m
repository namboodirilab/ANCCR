function [rpetimeline,valuetimeline,eventtimeline,inhibitiontimeline] =...
    simulatemicrostimulus(eventlog,rewardstate,nmicrostimulus,...
    statesize,alpha,gamma,lambda,sigma,d,inhibitionlog,maxstatelength)
% SIMULATEMICROSTIMULUS: simulate microstimulus TDRL model. 
%   eventlog: [eventindex, eventtime]
%   rewardstate: eventindex of rewards
%   nmicrostimulus: number of microstimuli each event elicits
%   statesize: size of state
%   alpha: learning rate
%   gamma: temporal discounting parameter
%   lambda: eligibility trace parameter
%   sigma: width of Gaussian function 
%   d: decay parameter
%   inhibitionlog: [start time, stop time, inhibition level (rpe)]
%           this allow to substitute RPE to artifically given inhibition 
%           level to mimic optogenetic inhibition of DA 
%   maxstatelength: truncate states at maxstatelength from each cue onset

if nargin<11
    maxstatelength = nan;    
end
nstimuli = length(unique(eventlog(:,1)));
sessionendtime = eventlog(end,2)+5;
nstate = floor(sessionendtime/statesize);

%% generate eventtimeline
% eventtimeline: [state index, reward magnitude]
% if reward happens multiple times within a single state, reward magnitude
% increases linearly
eventtimeline = zeros(nstate,2);
for is = 1:nstimuli
   eventidx = ceil(eventlog(eventlog(:,1)==is,2)/statesize);
   eventtimeline(unique(eventidx),1) = is;
   eventtimeline(unique(eventidx),2) = cellfun(@(x) sum(eventidx==x),num2cell(unique(eventidx)));
end

%% calculate inhibitiontimeline
% inhibitiontimeline shows level of inhibition in each time bin. If it is
% nan, there's no inhibition
inhibitiontimeline = nan(nstate,1);
for iinh = 1:size(inhibitionlog,1)
    ininh = ceil(inhibitionlog(iinh,1)/statesize):floor(inhibitionlog(iinh,2)/statesize);
    if size(inhibitionlog,2)==3
        inhibitiontimeline(ininh) = inhibitionlog(iinh,3);
    else
        inhstep = (inhibitionlog(iinh,4)-inhibitionlog(iinh,3))/(length(ininh)-1);
        inhibitiontimeline(ininh) = inhibitionlog(iinh,3)+[0:1:length(ininh)-1]*inhstep;
    end
end

%% calculate RPE 
[valuetimeline,rpetimeline] = deal(zeros(nstate,1));
%weight (w), eligibility trace (e), microstimulus level (x), memory trace height (y)
[w,e,x] = deal(zeros(nmicrostimulus,nstimuli));
y = deal(zeros(1,nstimuli));

% for first time bin 
stimulusidx = cell(nstimuli,1);
stimulusnum = nan(nstimuli,1);
for is = 1:nstimuli
   stimulusidx{is} = find(eventtimeline(:,1)==is); 
   stimulusnum(is) = sum(eventtimeline(:,1)==is);
end
tempstimuli = zeros(1,nstimuli);

rewarded = ismember(eventtimeline(1,1),rewardstate);
if rewarded
   reward = eventtimeline(1,2);
else
    reward = 0;
end
RPE = reward+gamma*valuetimeline(1);
if ~isnan(inhibitiontimeline(1))
    RPE = inhibitiontimeline(1);
end
rpetimeline(1) = RPE;

for is = 1:nstimuli
    w(:,is) = w(:,is)+alpha*rpetimeline(1)*e(:,is);
    e(:,is) = gamma*lambda*e(:,is)+x(:,is);
end

% iterate calculation for rest time bins
for i = 2:nstate
    skip = 0;

    % activate a set of micsotimuli when cs or us is delivered
    vtemp = 0;
    for is = 1:nstimuli
        if tempstimuli(is)+1<=stimulusnum(is)
            if i>=stimulusidx{is}(tempstimuli(is)+1)
                tempstimuli(is) = tempstimuli(is)+1;
            end
        end
        if tempstimuli(is)>0
            if ~isnan(maxstatelength)
                % if time longer than maxstatelength has passed since last
                % event, truncate state by skipping calculation
                if (i-stimulusidx{is}(tempstimuli(is)))*statesize>maxstatelength
                    valuetimeline(i) = valuetimeline(i-1);
                    skip = 1;
                    continue;
                end
            end
           y(is) = d^(i-stimulusidx{is}(tempstimuli(is)));
           x(:,is) = (1/sqrt(2*pi))*exp(-((y(is)-[1:nmicrostimulus]'/nmicrostimulus).^2)/(2*sigma^2))*y(is);
        end
        vtemp = vtemp + w(:,is)'*x(:,is);
    end
    if skip==1
        continue;
    end
    valuetimeline(i) = vtemp;

    % reward: reward magnitude of next time bin
    rewarded = ismember(eventtimeline(i,1),rewardstate);
    if rewarded
        reward = eventtimeline(i,2);
    else
        reward = 0;
    end
    
    RPE = reward+gamma*valuetimeline(i)-valuetimeline(i-1);
    if ~isnan(inhibitiontimeline(i))
        RPE = inhibitiontimeline(i);
    end
    rpetimeline(i) = RPE;
    
    % update weights and eligibility traces for cs and us
    for is = 1:nstimuli
       w(:,is) = w(:,is)+alpha*rpetimeline(i)*e(:,is);
       e(:,is) = gamma*lambda*e(:,is)+x(:,is);
    end
end
end




