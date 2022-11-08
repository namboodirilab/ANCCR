%Conditioned inhibition
clc; clearvars; close all;

%% Task setup
% Task parameter
cueITI = 30;
max_ITI = cueITI*3;
cueinhibperiod = 10;
cue_rew_delay = 3;
post_rew_delay = 3;
inter_cue_delay = 0;
numcue = 1000;
numtest = 3;
numrewards = 2000;
rewmag = 1;
rewITI = 5;

% Event labels
rew_label = 3;
cueinhibitor_label = 1;
cueinhibited_label = 2;

% Model parameters
samplingperiod = 0.2;
samplinginterval = samplingperiod;
alpha = 0.02;
alpha_r = alpha*10;
w = 0.5;
k = 1;
Tratio = 1.2;
minimumrate = 10^(-3);
threshold = 0.8;
maximumjitter = 0.1;
beta = [0; 0; 1];
Rtrue = [0; 0; 1];
n_stimuli = length(beta);
avg_window = 10;
exact_mean_or_not = 0;

nIter = 20;
for iiter = 1:nIter
    %% Simulate events
    % Negative association of CS2-reward
    eventlog_bkgd = simulateBackgroundRewards(numrewards, rewITI, rew_label, rewmag);
    eventlog_cue = simulateBackgroundRewards(numcue, cueITI, cueinhibitor_label, 0);
    eventlog_cue(:,2) = eventlog_cue(:,2);
    for i = 1:numcue % remove rewards that happen within [0 cueinhibperiod] from each cue
        rmidxs = find(eventlog_bkgd(:,2)>=eventlog_cue(i,2) & eventlog_bkgd(:,2)<eventlog_cue(i,2)+cueinhibperiod);
        eventlog_bkgd(rmidxs, :) = [];
    end
    eventlog_neg = [eventlog_bkgd;eventlog_cue];
    eventlog_neg = sortrows(eventlog_neg,2);
    eventlog_neg(find(eventlog_neg(:,1)==rew_label,1,'last')+1:end,:) = [];
    T = repmat(rewITI,size(eventlog_neg,1),1);
    
    % Positive association of CS1-reward
    eventlog_cs1 = simulateEvents(numcue, cueinhibited_label, rew_label, 1,...
        nan, cueITI, max_ITI, cue_rew_delay, rewmag, post_rew_delay);
    eventlog_cs1(:,2) = eventlog_cs1(:,2) + eventlog_neg(end,2);
    eventlog_pos = [eventlog_pos; eventlog_cs1];
    
    % Co-presentation of CS1 and CS2
    eventlog_copresent = simulateEventChain(numtest, 2, nan, cueITI, max_ITI, inter_cue_delay,...
        cue_rew_delay, 0, post_rew_delay);
    eventlog_copresent(:,2) = eventlog_copresent(:,2) + eventlog_pos(end,2);
    eventlog = [eventlog_pos; eventlog_copresent];
    T = [T; repmat(cueITI+post_rew_delay+cue_rew_delay,numcue*2+numtest*2,1)];
    
    %% Simulate ANCCR
    % Calculate model values
    [DA,ANCCR,PRC,SRC,NC,Rs] = calculateANCCR(eventlog, T*Tratio, alpha, k, ...
        samplingperiod, w,threshold,minimumrate,beta, alpha_r, maximumjitter, ...
        nan, nan);
    
    % Calculate values for C1 and C2
    ininhibitor = eventlog(:,1)==cueinhibitor_label; % C1
    value_inhibitor = squeeze(SRC(cueinhibitor_label,rew_label,ininhibitor).*Rs(cueinhibitor_label,rew_label,ininhibitor));
    
    ininhibited = eventlog(:,1)==cueinhibited_label; % C2
    value_inhibited = squeeze(SRC(cueinhibited_label,rew_label,ininhibited).*Rs(cueinhibited_label,rew_label,ininhibited));
    
    % Save values for C1 and C2, separate and paired
    value{1}(iiter,:) = value_inhibitor(end-numtest+1:end)+value_inhibited(end-numtest+1:end);
    value{2}(iiter,:) = value_inhibited(end-numtest*2+1:end-numtest);
    value{3}(iiter,:) = value_inhibitor(end-numtest*2+1:end-numtest);

    % Calculate action probabilities
    temperature  = 0.5;
    denom = exp(value{1}(iiter,:)/temperature) + exp(value{2}(iiter,:)/temperature) + exp(value{3}(iiter,:)/temperature);
    action_probs{1}(iiter,:) = exp(value{1}(iiter,:)/temperature) ./ denom;
    action_probs{2}(iiter,:) = exp(value{2}(iiter,:)/temperature) ./ denom;
    action_probs{3}(iiter,:) = exp(value{3}(iiter,:)/temperature) ./ denom;
end
%% Plot results
clr = {'k','r',;[0.6 0.6 0.6],[1 0.6 0.6]};

fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 2.5 3.5]);
axes('Position',axpt(5,5,2:5,1:4)); hold on;
for i = 2:3
    data_aps = mean(action_probs{i},2);
    bar(i-1,mean(data_aps),0.5,'FaceColor',clr{2,i-1});
    scatter(rand(nIter,1)*0.4-0.2+i-1,data_aps,5,clr{1,i-1},'filled');
    errorbar(i-1,mean(data_aps),std(data_aps)/sqrt(nIter),'k');
end

set(gca,'Box','off','TickDir','out','FontSize',8,'LineWidth',0.35,...
    'XTick',[1,2], 'XTickLabel',{'CS','CS + CI'},'XTickLabelRotation',45);
ylabel('p(lick|CS)')

savefig(fHandle,'conditionedinhibition.fig')
print(fHandle,'-depsc','-painters','conditionedinhibition.ai')
save('conditionedinhibition.mat', 'data_val', 'value', 'data_aps', ...
    'action_probs')

% T-test
CS = mean(action_probs{2},2);
CS_CI = mean(action_probs{3},2);
[significant_binary, p_value, ~, stats] = ttest(CS, CS_CI);
close all
