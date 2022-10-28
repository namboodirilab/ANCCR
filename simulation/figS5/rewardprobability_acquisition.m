clearvars; clc; close all;
rng(2);

%% parameter set up
% task parameters
numcue = 3000;
cuerewdelay = 1;
postrewdelay = 0;
rew_probs = 0.1:0.1:1;
meanITI = 10;

% anccr model parameters
samplingperiod = 0.2;   
alpha_anccr = 0.02;        
alpha_r = 0.2;
w = 0.5;               
k = 1;                 
minimumrate = 10^(-3);
maximumjitter = 0.1;
beta = [0,1];
threshold = 0.6;
Tratio = 1.2;

% rpe model parameters - csc/microstimulus
alpha_rpe = 0.001;
gamma = 0.95;
lambda = 0;
statesize = 1;

acqthreshold = [0.08; 0.6];
nIter = 100;
[acqtrialnum,acqrewardnum] = deal(nan(2,nIter,length(rew_probs)));
%%
for iIter = 1:nIter
    iIter
    for iP = 1:length(rew_probs)
        IRI = (meanITI+cuerewdelay+postrewdelay)/rew_probs(iP);
        eventlog_rr = simulateBackgroundRewards(500,IRI,2,1,1);
        eventlog = simulateEvents(numcue,1,2,...
            1,nan,meanITI,meanITI*3,cuerewdelay,rew_probs(iP),postrewdelay);
        for imdl = 1:2
            switch imdl
                case 1
                    [DA,~,eventtimeline] = simulateCSC(eventlog,2,statesize,alpha_rpe,gamma,lambda,[]);
                    incue = find(eventtimeline(:,1)==1);
                    inrw = find(eventtimeline(:,1)==2);
                case 2
                    eventlog(:,2) = eventlog(:,2)+eventlog_rr(end,2);
                    eventlog = [eventlog_rr;eventlog];
                    [DA,~,prc,src,nc,r] = calculateANCCR(eventlog, IRI*Tratio, alpha_anccr, k,...
                        samplingperiod,w,threshold,minimumrate,beta,alpha_r,maximumjitter,nan,nan);
                    incue = find(eventlog(:,1)==1,numcue,'last');
                    inrw = find(eventlog(:,1)==2);
                    inrw(1:500) = [];                    
            end
            cuersp = DA(incue);
            acqtrialnum(imdl,iIter,iP) = find(movmean(cuersp,30)>acqthreshold(imdl),1,'first');
            acqrewardnum(imdl,iIter,iP) = sum(inrw<incue(acqtrialnum(imdl,iIter,iP)));
        end
    end
end

%%
cd('D:\OneDrive - University of California, San Francisco\figures\manuscript\dopamine_contingency\revision\data');
save('rewardprobability_acquisition.mat','acqtrialnum','acqrewardnum','rew_probs');

%%

fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 7 4]);

for imdl = 1:2
    subplot(1,2,imdl)
    hold on;
     plot(log10(1./rew_probs),log10(squeeze(acqtrialnum(imdl,:,:)))','Color',[0.6 0.6 0.6]);
     plot(log10(1./rew_probs),log10(squeeze(acqrewardnum(imdl,:,:)))','Color',[1 0.6 0.6]);
     
     errorbar(log10(1./rew_probs),nanmean(log10(squeeze(acqtrialnum(imdl,:,:)))),...
         nanstd(log10(squeeze(acqtrialnum(imdl,:,:))))/sqrt(nIter),'color','k','CapSize',3);
     errorbar(log10(1./rew_probs),nanmean(log10(squeeze(acqrewardnum(imdl,:,:)))),...
         nanstd(log10(squeeze(acqrewardnum(imdl,:,:))))/sqrt(nIter),'color','r','CapSize',3);
     set(gca,'Box','off','TickDir','out','FontSize',8,'XLim',[-0.1 1.1],'XTick',log10(1./[1 0.5 0.3 0.1]),...
         'XTickLabel',{'1/1','2/1','3/1','10/1'},'YLim',[1 3.5],...
         'YTick',log10([20 50 100 500 1000 2000]),...
         'YTickLabel',{'20','50','100','500','1000','2000'},'YLim',log10([50 3000]));
     if imdl==1
         ylabel({'Predicted number required';'to acquisition (log scale)'});
         title({'RPE';'(model 1)'});
     else
         title({'ANCCR';'(model 2)'});
     end
end

%%
dir = 'D:\OneDrive - University of California, San Francisco\figures\manuscript\dopamine_contingency\revision';
print(fHandle,'-depsc','-painters',[dir,'\rewardprobability_acquisition.ai']);
