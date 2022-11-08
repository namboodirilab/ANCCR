% simulate different length of cue-reward delay with scaled or fixed ITI

clearvars; clc; close all;
rng(2);

%% parameter set up
% task parameters
numcue = 1000;
cuerewdelay = 3:10;
postrewdelay = 0;
rew_probs = 1;
Tratio = 1.2;

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

% rpe model parameters - csc
alpha_rpe = 0.025;
gamma = 0.95;
lambda = 0;
statesize = 1;

% threshold defining acqusition trial
% acquisition trial was defined as the first trial when DA cue response 
% exceeds this threshold
acqthreshold = [0.08; 0.6];
nIter = 100;

acqtrialnum = nan(2,nIter,length(cuerewdelay),2);
%%
for iIter = 1:nIter
    iIter
    for iD = 1:length(cuerewdelay)
        for iI = 1:2
            if iI==1
                % fixed ITI: ITI was set to 30s regardless of length of cue-reward delay
                meanITI = cuerewdelay(1)*10; 
            else
                % scaled ITI: ITI was set to 10 times of cue-reward delay
                meanITI = cuerewdelay(iD)*10;
            end
            % generate eventlog
            IRI = meanITI+cuerewdelay(iD)+postrewdelay;
            eventlog = simulateEvents(numcue,1,2,...
                1,nan,meanITI,meanITI*3,cuerewdelay(iD),rew_probs,postrewdelay);

            for imdl = 1:2
                switch imdl
                    case 1
                        [DA,~,eventtimeline] = simulateCSC(eventlog,2,statesize,alpha_rpe,gamma,lambda,[]);
                        incue = eventtimeline(:,1)==1;
                    case 2
                        DA = calculateANCCR(eventlog, IRI*Tratio, alpha_anccr, k,...
                            samplingperiod,w,threshold,minimumrate,beta,alpha_r,maximumjitter,nan,nan);
                        incue = eventlog(:,1)==1;
                end

                % calculate number of trials until acquisition
                cuersp = DA(incue);
                acqtrialnum(imdl,iIter,iD,iI) = find(movmean(cuersp,30)>acqthreshold(imdl),1,'first');
            end
        end
    end
end

%% save data
cd('D:\OneDrive - University of California, San Francisco\figures\manuscript\dopamine_contingency\revision\data');
save('scaledITI_acquisition.mat','acqtrialnum','cuerewdelay','acqthreshold');

%% FigS5C
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 7 4]);

for imdl = 1:2
    subplot(1,2,imdl)
    hold on;
    plot(cuerewdelay,squeeze(acqtrialnum(imdl,:,:,1)),'Color',[0.6 0.6 0.6]);
    plot(cuerewdelay,squeeze(acqtrialnum(imdl,:,:,2)),'Color',[1 0.6 0.6]);
    
    errorbar(cuerewdelay,mean(squeeze(acqtrialnum(imdl,:,:,1))),...
        std(squeeze(acqtrialnum(imdl,:,:,1)))/sqrt(nIter),'color','k','CapSize',3);
    errorbar(cuerewdelay,mean(squeeze(acqtrialnum(imdl,:,:,2))),...
        std(squeeze(acqtrialnum(imdl,:,:,2)))/sqrt(nIter),'color','r','CapSize',3);
    set(gca,'Box','off','TickDir','out','FontSize',8,'XLim',[2 11],'XTick',3:3:9,...
        'YLim',[0 400],'YTick',0:200:400);
     if imdl==1
         ylabel({'Predicted trials';'to acquisition'});
         title({'RPE';'(model 1)'});
     else
         title({'ANCCR';'(model 2)'});
     end
end
