function [rpetimeline,valuetimeline,eventtimeline,statetimeline,inhibitiontimeline] =...
    simulateCSC(eventlog,rewardstate,statesize,...
    alpha,gamma,lambda,inhibitionlog,withITIstates,cuetoITIdelay,maxstatelength)

if nargin<8
   withITIstates = 1; 
end
if nargin<10
    maxstatelength = nan;
end
%inhibitionlog: start time, stop time, inhibition level (rpe)

nstimuli = length(unique(eventlog(:,1)));
states = cell(nstimuli,1);
sessionendtime = eventlog(end,2)+5;
nstate = floor(sessionendtime/statesize);
% 1st column: eventidx
% 2nd column: number of events in each state; to account for when multiple rewards happen in a single state
eventtimeline = zeros(nstate,1); 
statetonext = [diff(eventlog(:,2));nan];
statespacesize = 0;

for is = 1:nstimuli
    if is==rewardstate
        numstates = ceil(max([statetonext(eventlog(:,1)==is);eventlog(1,2)])/statesize);
    else
        numstates = ceil(max(statetonext(eventlog(:,1)==is))/statesize);
    end
    if ~isnan(maxstatelength)
    numstates = min([numstates,ceil(maxstatelength/statesize)]);
    end
    states{is} = [1:numstates]+statespacesize;
    statespacesize = statespacesize+numstates;
    eventidx = ceil(eventlog(eventlog(:,1)==is,2)/statesize);
    eventtimeline(unique(eventidx),1) = is;
    eventtimeline(unique(eventidx),2) = cellfun(@(x) sum(eventidx==x),num2cell(unique(eventidx)));
end

if withITIstates==0
   ITIstates = [states{rewardstate}];
   for is = unique(eventlog(:,1))'
      if is==rewardstate
          continue;
      end
      ITIstates = [ITIstates, states{is}(ceil(cuetoITIdelay(is)/statesize)+1:end)];
   end
else
    ITIstates = [];
end

statetimeline = NaN(nstate, 1);
temp1 = ceil(eventlog(1,2)/statesize)-1;
ns = min([length(states{rewardstate}),temp1]);
if temp1>length(states{rewardstate})
    statetimeline(1:temp1) = [states{rewardstate}(1:ns),repmat(states{rewardstate}(end),1,temp1-ns)];
else
    statetimeline(1:temp1) = states{rewardstate}(1:temp1);
end
temp = temp1+1;
nevent = size(eventlog,1);
for ie = 2:nevent
    temp1 = ceil(eventlog(ie,2)/statesize)-1;
    ns = min([length(states{eventlog(ie,1)}),temp1-temp+1]);
    statetimeline(temp:temp1) = [states{eventlog(ie-1,1)}(1:ns),repmat(states{eventlog(ie-1,1)}(end),1,temp1-temp+1-ns)];
    temp = temp1+1;
end
temp1 = nstate-temp+1;
ns = min([temp1,length(states{eventlog(nevent,1)})]);
statetimeline(temp:end) = [states{eventlog(nevent,1)}(1:ns),repmat(states{eventlog(nevent,1)}(end),1,temp1-ns)];

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

[valuetimeline,rpetimeline] = deal(zeros(nstate,1));
[valuesforstate,eforstate] = deal(zeros(statespacesize,1));

rewarded = ismember(eventtimeline(1,1),rewardstate);
if rewarded
    reward = eventtimeline(1,2);
else
    reward = 0;
end
RPE = reward + gamma*valuesforstate(statetimeline(1));
rpetimeline(1) = RPE;
eforstate = gamma*lambda*eforstate;
eforstate(statetimeline(1)) = eforstate(statetimeline(1))+1;

for i = 1:nstate-1
    if statetimeline(i)==statetimeline(i+1) && ismember(statetimeline(i),cellfun(@(x) x(end),states)) && ~isnan(maxstatelength)
        valuetimeline(i) = valuetimeline(i-1);
        continue
    end
    
    if ismember(statetimeline(i),ITIstates) & withITIstates==0
        v = 0;
    else
        v = valuesforstate(statetimeline(i));
    end
    
    if ismember(statetimeline(i+1),ITIstates) & withITIstates==0
        v1 = 0;
    else
        v1 = valuesforstate(statetimeline(i+1));
    end
    
    rewarded = ismember(eventtimeline(i+1,1),rewardstate);
    if rewarded
        reward = eventtimeline(i+1,2);
    else
        reward = 0;
    end
    
    RPE = reward + gamma*v1 - v;
    if ~isnan(inhibitiontimeline(i+1))
        RPE = inhibitiontimeline(i+1);
    end
    rpetimeline(i+1) = RPE;

    if ismember(statetimeline(i),ITIstates) & withITIstates==0
        continue;
    end
    eforstate = gamma*lambda*eforstate;
    eforstate(statetimeline(i)) = eforstate(statetimeline(i))+1;

    valuesforstate = valuesforstate + alpha*RPE*eforstate;
    valuetimeline(i) = valuesforstate(statetimeline(i));
end
end