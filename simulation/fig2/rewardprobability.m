%DA signals +/- prediction errors to US and CS, 
% in proportion to reward magnitude and/or probability of reward.
clearvars; clc; 

rng(2)
%% Task setup
% Task parameters
meanITI = 100;
maxITI = meanITI*3;
cuerewdelay = 0.5;
postrewdelay = 1;
numcue = 3000;

% Model parameters
samplingperiod = 0.2;
alpha = 0.02;
alpha_r = 0.2;
w = 0.6;
k = 0.01;
Tratio = 1.2;
minimumrate = 10^(-3);
threshold = 0.6;
maximumjitter = 0.1;
beta = [0; 1; 0];

nIter = 20;
rew_probs = [0.3; 0.6; 0.9];
IRI = (meanITI + cuerewdelay) ./ rew_probs;

probability_results = nan(nIter, length(rew_probs), length(beta));

%% Run simulation
% Simulate for cue reward pairs with three different reward
% probabilities
DA = [];
for j = 1:length(rew_probs)
    for iIter = 1:nIter
        % Simulate cue/reward delivery
        [eventlog] = simulateEvents(numcue, 1, 2, 1, 3,  ...
            meanITI, maxITI, cuerewdelay, rew_probs(j), postrewdelay);
        
        om_resp = zeros(sum(eventlog(:,1)==3),1);
        omidx = find(eventlog(:,1)==3);
        % Recognition of omission state is dependent on the probability of reward
        % if the reward probability is 10%, you would expect omission state
        % after cue only in 10% probability
        outom = sort(randsample(length(omidx),round(length(omidx)*(1-rew_probs(j)))));
        eventlog(omidx(outom),:) = [];
         
        % Calculate model values
        [DA,ANCCR,PRC,SRC,NC,Rs] = calculateANCCR(eventlog, IRI(2)*Tratio, alpha, k, ...
            samplingperiod ,w,threshold,minimumrate,beta,alpha_r,maximumjitter,nan,[3,2]);
        
        % Save cue, reward, omission indices
        cue_resp = DA(eventlog(:,1) == 1);
        rew_resp = DA(eventlog(:,1) == 2);
        om_resp(~ismember(1:length(omidx),outom)) = DA(eventlog(:,1) == 3);
        probability_results(iIter, j, 1) = mean(cue_resp(end-100:end));
        probability_results(iIter, j, 2) = mean(rew_resp(end-100:end));
        probability_results(iIter, j, 3) = mean(om_resp(end-100:end));
   end
end

%% Plotting DA response

close all
dir = 'D:\OneDrive - University of California, San Francisco\figures\manuscript\dopamine_contingency\revision';

% Plot reward probability results
clr_light = {[0.6 0.6 0.6],[0.6 0.6 1],[1 0.6 0.6]};
clr = {[0 0 0],[0 0 1],[1 0 0]};
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 2.5 3.5]);
axes('Position',axpt(5,5,2:5,1:4)) 
hold on;
for i = 1:3
plot(1:3,squeeze(probability_results(:,:,i)),'Color',clr_light{i},'LineWidth',0.35);
errorbar(1:3,mean(squeeze(probability_results(:,:,i)),1),std(squeeze(probability_results(:,:,i)),[],1)/sqrt(nIter),...
    'Color',clr{i},'LineWidth',0.5);
end
plot([0 4],[0 0],'k:','LineWidth',0.35);
set(gca,'Box','off','TickDir','out','FontSize',8,'LineWidth',0.35,...
    'XTick',1:3,'XTickLabel',rew_probs*100,...
    'YTick',-0.5:0.5:1,'YLim',[-0.3 1.3],'XLim',[0.5 3.5])
xlabel('Reward probability (%)')
ylabel('Predicted DA response')

print(fHandle,'-depsc','-painters',[dir,'\rew_prob.ai']);

