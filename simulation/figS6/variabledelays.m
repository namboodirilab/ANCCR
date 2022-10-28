clearvars; clc; close all;
rng(2);

%% parameter set up
% task parameters
numcue = 1000;
cuerewdelay = [1 10 100];
meanITI = 50;
rew_probs = 1;
IRI = meanITI/3;

% anccr model parameters
samplingperiod = 0.2;   
alpha_anccr = 0.02;        
alpha_r = 0.2;
w = 0.5;               
k = 1;                 
minimumrate = 10^(-3);
maximumjitter = 0.1;
beta = [0,0,0,1,1,1];
threshold = 0.6;
Tratio = 5;
exact_mean_or_not = 1;

nIter = 100;
%%
avecuersp = nan(nIter,3,2);
for iIter = 1:nIter
    iIter
    for iI = 1:2
        if iI==1
            eventlog = simulateEventsTrialLess(repmat(numcue,1,3),1:3,4:6,[1,1,1],...
                [nan,nan,nan],repmat(meanITI,1,3),repmat(meanITI*3,1,3),[0,0,0],cuerewdelay,[1,1,1]);
        else
            eventlog = simulateBackgroundRewards(repmat(numcue,1,6),....
                repmat(meanITI,1,6),1:6,[zeros(1,3),ones(1,3)],1);
        end
                
        [DA,ANCCR,~,~,NC] = calculateANCCR(eventlog, IRI*Tratio, alpha_anccr, k,...
            samplingperiod,w,threshold,minimumrate,beta,alpha_r,maximumjitter,nan,nan,exact_mean_or_not);
        for icue = 1:length(cuerewdelay)
            incue = eventlog(:,1)==icue;
            cuersp = DA(incue);
            avecuersp(iIter,icue,iI) = mean(cuersp(end-99:end));
        end
    end
end

%%
cd('D:\OneDrive - University of California, San Francisco\figures\manuscript\dopamine_contingency\revision\data');
save('variabledelays.mat','avecuersp','cuerewdelay','meanITI');

%%
clr_light = {[1 0.8 0.8],[0.8 0.8 0.8]};
clr = {[1 0 0],[0 0 0]};
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 3 4]);
hold on;
plot([0 4],[0 0],'k:','LineWidth',0.35);
plot(1:3,avecuersp(:,:,1),'Color',clr_light{1},'LineWidth',0.35);
plot(1:3,avecuersp(:,:,2),'Color',clr_light{2},'LineWidth',0.35);
errorbar(1:3,mean(avecuersp(:,:,1)),std(avecuersp(:,:,1))/sqrt(nIter),'Color',clr{1},'LineWidth',0.5);
errorbar(1:3,mean(avecuersp(:,:,2)),std(avecuersp(:,:,2))/sqrt(nIter),'Color',clr{2},'LineWidth',0.5);
set(gca,'XTick',1:3,'XTickLabel',{'C1';'C2';'C3'},'XTickLabelRotation',45,...
    'YTick',0:0.5:1,'Box','off','TickDir','out',...
    'FontSize',8,'LineWidth',0.35,'XLim',[0.5 3.5],'YLim',[-0.2 1.1]);
ylabel('Predicted cue response');

%%
dir = 'D:\OneDrive - University of California, San Francisco\figures\manuscript\dopamine_contingency\revision';
print(fHandle,'-depsc','-painters',[dir,'\variabledelays.ai']);
