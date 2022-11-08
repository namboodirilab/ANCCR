% simulate trial-less task with various cue-reward delays using CSC and
% ANCCR models
clearvars; clc; close all;
rng(2);

%% parameter set up
% task parameters
numcue = 2000;
meanITI = 12;
cuerewdelay = [0.05,0.1,0.2,0.5,1,2,5]*meanITI;
rew_probs = 1;

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
Tratio = 5;
exact_mean_or_not = 1;

% rpe model parameters - csc/microstimulus
alpha_rpe = 0.05;
gamma = 0.95;
lambda = 0;
statesize = 0.2;

nIter = 100;
%%
avecuersp = nan(2,length(cuerewdelay)+1,nIter);
for iIter = 1:nIter
    iIter
    for iD = 1:length(cuerewdelay)+1
        
        if iD<length(cuerewdelay)+1
            % cue-reward association w/ variable cue-reward delays
            eventlog = simulateBackgroundRewards(numcue,....
                repmat(meanITI,1,2),1,0,1);
            rwtimes = eventlog(:,2)+cuerewdelay(iD);
            eventlog = [eventlog;[ones(numcue,1)*2,rwtimes,ones(numcue,1)]];
            eventlog = sortrows(eventlog,2);
        else
            % no cue-reward association control
            eventlog = simulateBackgroundRewards(repmat(numcue,1,2),....
                repmat(meanITI,1,2),[1,2],[0,1],1);
        end
        for imdl = 1:2
            switch imdl
                case 1
                    [DA,~,eventtimeline] = simulateCSC(eventlog,2,statesize,alpha_rpe,gamma,lambda,[]);
                    incue = eventtimeline(:,1)==1;
                case 2
                    [DA,ANCCR,~,~,NC] = calculateANCCR(eventlog, meanITI*Tratio, alpha_anccr, k,...
                        samplingperiod,w,threshold,minimumrate,beta,alpha_r,maximumjitter,nan,nan,exact_mean_or_not);
                    incue = eventlog(:,1)==1;
            end
            % calculate averaged cue response over last 500 cues
            cuersp = DA(incue);
            avecuersp(imdl,iD,iIter) = mean(cuersp(end-499:end));
        end
    end
end

%% save data
cd('D:\OneDrive - University of California, San Francisco\figures\manuscript\dopamine_contingency\revision\data');
save('timescale.mat','avecuersp','cuerewdelay','meanITI');

%% Fig6SB
ct = cbrewer('qual','Dark2',3);
clr = [0 0 0; 0.6 0.6 0.6];

fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 6 4.5]);
for imdl = 1:2
    subplot(1,2,imdl)
    hold on;
    data = squeeze(avecuersp(imdl,1:end-1,:));
    data_ctrl = sort(squeeze(avecuersp(imdl,end,:)));

    plot(cuerewdelay/meanITI,data,'Color',[0.6 0.6 0.6],'LineWidth',0.35);
    errorbar(cuerewdelay/meanITI,nanmean(data,2),nanstd(data,[],2)/sqrt(nIter),'Color','k','LineWidth',0.5);
    plot([-1 6],repmat(nanmean(data_ctrl),1,2),'k--');
    
    set(gca,'XTick',[0 0.5 1 2 5],'XTickLabelRotation',45,...
        'YTick',0:0.5:1,'Box','off','TickDir','out',...
        'FontSize',8,'LineWidth',0.35,'XLim',[-0.5 5.5],'YLim',[-0.1 1]);
    if imdl==1
        ylabel('Predicted cue response');
        title({'RPE';'(model 1)'});
    else
        title({'contingency';'(model 2)'});
    end
    
end
