function [DA,ANCCR,PRC,SRC,NC,Rs,Delta,Mij,Mi] =...
    calculateANCCR(eventlog, T, alpha, k,samplinginterval,w,threshold,...
    minimumrate,beta,alpha_r,maximumjitter,optolog,omidx,exact_mean_or_not,nevent_for_edge)
if alpha_r>1
    alpha_r = 1;
end

if nargin<=11 | isnan(optolog)
    optolog = zeros(size(eventlog,1),2);
end
if nargin<=12 | isnan(omidx)
    % First entry is omission index, second entry is corresponding reward index
    omidx = [nan,nan];
end
if nargin<=13
    % if exact_mean_or_not=1, calculate exact mean for Mij instead using
    % alpha
    exact_mean_or_not = 0;
end

if nargin<=14
    % if nevent_for_edge>0, use averaged NC for last nevent to calculate
    % edge
    nevent_for_edge = 0; 
end

% omtrue: whether the omission state will be used or not in the calculaton of ANCCR
omtrue = false(size(omidx,1),1);
uniquetime = unique(eventlog(:,2));

%% if more than one event happens at the same time, assume random perceptual delay between them
for jt = 1:length(uniquetime)
    if sum(eventlog(:,2)==uniquetime(jt))==1
        continue;
    end
    idx = find(eventlog(:,2)==uniquetime(jt));
    eventlog(idx(2:end),2) = eventlog(idx(2:end),2)+randn(length(idx)-1,1)*maximumjitter;
end
eventlog = sortrows(eventlog,2);
ntime = size(eventlog,1);

%%
nstimuli = length(unique(eventlog(:,1)));
samplingtime = 0:samplinginterval:eventlog(end,2);

% if T is a vector, use T(jt) for the calculation at time jt. othersiwse,
% use fixed T
if length(T)==1
    T = repmat(T,size(eventlog,1),1);
end
gamma = exp(-1./T);

[Eij,Ei,Mi,Delta] = deal(zeros(nstimuli,ntime));
[Mij,PRC,SRC,NC,ANCCR,Rs] = deal(zeros(nstimuli,nstimuli,ntime));
R = zeros(nstimuli,nstimuli);
numevents = zeros(nstimuli,1);
DA = zeros(ntime,1);

beta = beta(unique(eventlog(:,1)));
Imct = beta(:)>threshold;
nextt = 1;
numsampling = 0;
%%
for jt = 1:ntime
    skip = false;
    je = eventlog(jt,1);

    if ismember(je,omidx(:,1))
        if ~omtrue(omidx(:,1)==je)
            Delta(:,jt) = Delta(:,jt-1);
            Eij(:,jt) = Eij(:,jt-1);
            Mij(:,:,jt) = Mij(:,:,jt-1);
            PRC(:,:,jt) = PRC(:,:,jt-1);
            SRC(:,:,jt) = SRC(:,:,jt-1);
            NC(:,:,jt) = NC(:,:,jt-1);
            skip = true;
        end
    end
    
    if ~skip
        numevents(je) = numevents(je)+1;
        if exact_mean_or_not == 0            
            if ~isstruct(alpha)
                alphat = alpha;
            else
                % if alpha is structure, alpha exponentially decreases from
                % alpha.init to alpha.min w/ alpha.exponent decrease constant
                alphat = exp(-alpha.exponent*(jt-0))*(alpha.init-alpha.min)+alpha.min;
            end
        else
            alphat = 1/numevents(je);
        end
    
        if jt>1
            Delta(:,jt) = Delta(:,jt-1)*gamma(jt)^(eventlog(jt,2)-eventlog(jt-1,2));
            Eij(:,jt) = Eij(:,jt-1)*gamma(jt)^(eventlog(jt,2)-eventlog(jt-1,2));
            Mij(:,:,jt) = Mij(:,:,jt-1);
            ANCCR(~ismember(1:nstimuli,je),:,jt) = ANCCR(~ismember(1:nstimuli,je),:,jt-1);
        end
        % Indicator for whether event has recently happened
        % Delta resets to one at every instance of event w/o cumulative sum
        Delta(je,jt) = 1;
        Eij(je,jt) = Eij(je,jt)+1;
        Mij(:,je,jt) = Mij(:,je,jt)+alphat*(Eij(:,jt)-Mij(:,je,jt)).*Imct(je);
        
        PRC(:,:,jt) = Mij(:,:,jt)-repmat(Mi(:,jt),1,nstimuli);
        SRC(:,:,jt) = PRC(:,:,jt).*repmat(Mi(:,jt)',nstimuli,1)./repmat(Mi(:,jt),1,nstimuli);
        belowminrate = Mi(:,jt)/T(jt)<minimumrate;
        SRC(belowminrate,:,jt) = 0;
        
        % something to make sure only calculating contingency and R after experiencing
        % first Y; this part can be improved later
        PRC(numevents==0,:,jt) = 0;
        PRC(:,numevents==0,jt) = 0;
        SRC(numevents==0,:,jt) = 0;
        SRC(:,numevents==0,jt) = 0;
        R(:,numevents==0) = 0;
        R(numevents==0,:) = 0;
        
        NC(:,:,jt) = w*SRC(:,:,jt)+(1-w)*PRC(:,:,jt);
        
        % Indicator for whether an event is associated with another event
         Iedge = mean(NC(:,je,max([1,jt-nevent_for_edge]:jt)),3)>threshold;
        Iedge(je) = false;
        
        % once the cause of reward state is revealed, omission state of that
        % reward state can be used for calculation of ANCCR. Before that,
        % omission state is ignored
        if ismember(je,omidx(:,2)) && sum(Iedge)>0
            omtrue(omidx(:,2)==je) = omtrue(omidx(:,2)==je) | true;
        end
        
        % calculate ANCCR for every event
        % Rjj is externally driven; the magnitude of stimulus an animal just experienced
        R(je,je) = eventlog(jt,3);
            
        for ke = 1:nstimuli
            Iedge_ke = mean(NC(:,ke,max([1,jt-nevent_for_edge]:jt)),3)>threshold;
            Iedge_ke(ke) = false;
            
            ANCCR(ke,:,jt) = NC(ke,:,jt).*R(ke,:)-...
                    sum(ANCCR(:,:,jt).*Delta(:,jt).*repmat(Iedge_ke,1,nstimuli));
        end


        if ~(optolog(jt,1) == 1) % If target is not inhibited, normal DA
            DA(jt) = sum(ANCCR(je,:,jt).*Imct');
        else % If target is inhibited, replace DA
            DA(jt) = optolog(jt,2);
        end

        if ismember(je,omidx(:,1))
            je_om = find(je==omidx(:,1));
            % if the current state is omission of j, R(omission,j) = R(j,j)
            R(je,omidx(je_om,2)) = R(omidx(je_om,2),omidx(je_om,2));
            % omission state is an MCT
            Imct(je) = true;
        end
        % This must come after opto s.t. Imct is not formed before opto applied
        Imct(je) = Imct(je) | DA(jt)+beta(je)>threshold;

        Rs(:,:,jt) = R;
        if DA(jt)>=0
            R(:,je) = R(:,je)+alpha_r*(eventlog(jt, 3)-R(:,je));
        else
            if any(Iedge)
                R(Iedge,je) = R(Iedge,je) -...
                    alpha_r*R(Iedge,je).*((Delta(Iedge,jt)./numevents(Iedge)) ./ sum((Delta(Iedge,jt)./numevents(Iedge))));
            else
                R(:,je) = R(:,je);
            end
        end
        
    end
    
    if jt<ntime
        subsamplingtime = samplingtime(samplingtime>=eventlog(jt,2) & samplingtime<eventlog(jt+1,2));
              
        Ei(:,jt+1) = Ei(:,jt)*gamma(jt)^samplinginterval;
        if ~isempty(subsamplingtime)
            for jjt = nextt:jt
                if ismember(eventlog(jjt,1),omidx(:,1))
                    if ~omtrue(omidx(:,1)==eventlog(jjt,1))
                        continue
                    end
                end
                Ei(eventlog(jjt,1),jt+1) = Ei(eventlog(jjt,1),jt+1)+...
                    gamma(jt).^(subsamplingtime(1)-eventlog(jjt,2));
            end
            nextt = jt+1;
        end
        if exact_mean_or_not == 0
            if ~isstruct(alpha)
                alphat = alpha;
            else
                alphat = exp(-alpha.exponent*(jt-0))*(alpha.init-alpha.min)+alpha.min;
            end
        else
            alphat = 1/(numsampling+1);
        end
        
        Mi(:,jt+1) = Mi(:,jt)+k*alphat*(Ei(:,jt+1)-Mi(:,jt));
        for iit = 2:length(subsamplingtime)
            if exact_mean_or_not == 0
                if ~isstruct(alpha)
                    alphat = alpha;
                else
                    alphat = exp(-alpha.exponent*(jt-0))*(alpha.init-alpha.min)+alpha.min;
                end
            else
                alphat = 1/(numsampling+iit);
            end
        
            Ei(:,jt+1) = Ei(:,jt+1)*gamma(jt)^samplinginterval;
            Mi(:,jt+1) = Mi(:,jt+1)+k*alphat*(Ei(:,jt+1)-Mi(:,jt+1));
        end
        numsampling = numsampling+length(subsamplingtime); 
    end
end
end