clearvars; clc;
rng(2)

%% Parameter set up
% Task parameters
numcue = [500,100]; % 100% C-R conditioning; 50% C-R conditioning
cuerewdelay = 1;
postrewdelay = 0;
rew_probs = [1,0.5];
meanITI = 30;
Tratio = 1.2;
numrewards = 500; % bgd rewards before conditioning
IRI = meanITI+cuerewdelay+postrewdelay;

% Model parameters
samplingperiod = 0.2;   
alpha_anccr = 0.02;        
alpha_r = 0.2;
w = 0.6;               
k = 0.01;                 
minimumrate = 10^(-3);
maximumjitter = 0.1;
beta = [0,1,0];
threshold = 0.6;

% CSC model parameters
alpha_csc = 0.05;
gamma = 0.95;
lambda = 0;
statesize = 0.2;

%% Simulation

nIter = 20;
omave = nan(nIter,2);
[cueave,rewave] = deal(nan(nIter,2,2));
numtest = 5;
for iIter = 1:nIter
    % Simulate background rewards
    eventlog_bgd = simulateBackgroundRewards(numrewards,meanITI+cuerewdelay+postrewdelay,2,1,0);
    
    eventlog_pre = simulateEvents(numcue(1), 1, 2, 1, 3,  ...
        meanITI, meanITI*3, cuerewdelay, rew_probs(1), postrewdelay);
    eventlog_pre(:,2) = eventlog_pre(:,2)+eventlog_bgd(end,2);

    eventlog = simulateEvents(numcue(2), 1, 2, 1, 3, ...
        meanITI, meanITI*3, cuerewdelay, rew_probs(2), postrewdelay);
    eventlog(:,2) = eventlog(:,2)+eventlog_pre(end,2);

    eventlog = [eventlog_bgd;eventlog_pre;eventlog];

    for imdl = 1:2
        if imdl==1
            % Simulate CSC
            [DA,~,eventtimeline,statetimeline] = simulateCSC(eventlog,2,statesize,alpha_csc,gamma,lambda,[]);
            omrsp = DA(eventtimeline==3);        
        else
            % Simulate omission state
            [eventlog,cueoutcomelog] = simulateOmissionState(eventlog, 1, 3, rew_probs(2));
            cueoutcomelog = cueoutcomelog{1};
            % Calculate ANCCR model values
            DA = calculateANCCR(eventlog, IRI*Tratio, alpha_anccr, k, ...
                samplingperiod,w,threshold,minimumrate,beta,alpha_r,maximumjitter,nan,[3,2]);
            omrsp = zeros(sum(cueoutcomelog(:,2)<=0),1);
            omrsp(cueoutcomelog(cueoutcomelog(:,2)<=0,2)==0) = DA(eventlog(:,1)==3);
            eventtimeline = eventlog(:,1);
        end

        omave(iIter,imdl) = mean(omrsp);

        % Find cue indices and average response
        cueidx = find(eventtimeline(:,1)==1);
        cueave(iIter,imdl,1) = mean(DA(cueidx(1:numtest)));
        cueave(iIter,imdl,2) = mean(DA(cueidx([-numtest+1:0]+numcue(1))));

        % Find reward indices and average response
        rewidx = find(eventtimeline(:,1)==2);
        rewidx(rewidx<cueidx(1)) = []; % remove bgd rewards
        rewave(iIter,imdl,1) = mean(DA(rewidx(1:numtest)));
        rewave(iIter,imdl,2) = mean(DA(rewidx([-numtest+1:0]+numcue(1))));
    end
end

% Plot values
imdl = 2;
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 2.5 3.5]);
axes('Position',axpt(5,5,2:5,1:4)) 
hold on;
    plot(1:2,squeeze(cueave(:,imdl,:)),'Color',[0.6 0.6 0.6],'LineWidth',0.35);
    plot(1:3,squeeze([squeeze(rewave(:,imdl,:)),omave(:,imdl)]),'Color',[0.6 0.6 1],'LineWidth',0.35);
    errorbar(1:2,nanmean(squeeze(cueave(:,imdl,:))),nanstd(squeeze(cueave(:,imdl,:)))/sqrt(nIter),'k','LineWidth',0.5);
    errorbar(1:3,nanmean([squeeze(rewave(:,imdl,:)),omave(:,imdl)]),...
        nanstd([squeeze(rewave(:,imdl,:)),omave(:,imdl)])/sqrt(nIter),'b','LineWidth',0.5);
    plot([-0.5 3.5],[0 0],'k:');

    ylim([-0.3 1.3])
    xlim([0.5 3.5]);
    set(gca,'XTick',[1 2 3],'XTickLabel',{'early','late','omission'},'XTickLabelRotation',45,...
    'Box','off','TickDir','out','FontSize',8,'LineWidth',0.35,'YTick',-0.5:0.5:1.5)
    ylabel('Predicted DA response')

print(fHandle,'-depsc','-painters','D:\OneDrive - University of California, San Francisco\dopamine contingency\figures\figures for r3\omission.ai')

