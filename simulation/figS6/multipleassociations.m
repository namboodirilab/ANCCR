clearvars; clc; close all;
rng(2);

%% parameter set up
% task parameters
numcue = 1000;
cuerewdelay = [2,2,4,0,0,0,0]; %c1,2,3,5,6,7,8
cuecuedelay = 2;
meanITI = 30;
rew_probs = [1,1,1,0,0,0,0];
IRI = meanITI/3;

% anccr model parameters
samplingperiod = 0.2;   
alpha_anccr = 0.02;        
alpha_r = 0.2;
w = 0.5;               
k = 1;                 
minimumrate = 10^(-3);
maximumjitter = 0.1;
beta = [0,0,0,0,0,0,0,0,1,1,1];
threshold = 0.6;
Tratio = 1.2;
exact_mean_or_not = 1;

nIter = 100;
NCave = nan(nIter,length(beta),length(beta));
nave = 1000;
%%
avecuersp = nan(nIter,3);
for iIter = 1:nIter
    iIter
    eventlog = simulateEventsTrialLess(repmat(numcue,1,7),[1,2,3,5,6,7,8],[9,10,11,nan,nan,nan,nan],...
        [1,1,1,nan,nan,nan,nan],nan,meanITI,meanITI*3,0,cuerewdelay,rew_probs);
    cs3idx = eventlog(:,1)==3;
    eventlog_c4 = [ones(sum(cs3idx),1)*4,eventlog(cs3idx,2)+cuecuedelay,zeros(sum(cs3idx),1)];
    eventlog = [eventlog;eventlog_c4];
    eventlog = sortrows(eventlog,2);
    
    [DA,ANCCR,~,~,NC] = calculateANCCR(eventlog, IRI*Tratio, alpha_anccr, k,...
        samplingperiod,w,threshold,minimumrate,beta,alpha_r,maximumjitter,nan,nan,exact_mean_or_not);
    NCave(iIter,:,:) = squeeze(mean(NC(:,:,end-nave+1:end),3));
end

%%
cd('D:\OneDrive - University of California, San Francisco\figures\manuscript\dopamine_contingency\revision\data');
save('multipleassociations.mat','NCave');

%%

fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 5.8 4.2]);
imagesc(squeeze(mean(NCave)));
statelist = {'C1';'C2';'C3';'C4';'C5';'C6';'C7';'C8';'R1';'R2';'R3'};
set(gca,'XTick',1:11,'XTickLabel',statelist,'Box','off','TickDir','out','FontSize',8,...
    'YTick',1:11,'YTickLabel',statelist,'CLim',[-1 1],'XTickLabelRotation',90);
hc = colorbar;
set(hc,'XTick',-0.5:0.5:1,'Box','off','TickDir','out','XLim',[-0.5 1]);

%%
dir = 'D:\OneDrive - University of California, San Francisco\figures\manuscript\dopamine_contingency\revision';
print(fHandle,'-depsc','-painters',[dir,'\multipleassociations.ai']);
