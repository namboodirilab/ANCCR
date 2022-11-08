% simulate serialcue conditioning with ANCCR model
% dopamine response was inhibited during cs2 over learning

clearvars; close all; clc;
rng(2)

% task parameter
numcue = 400;
cuerewdelay = 1.5; % delay from cs2 to reward
cuecuedelay = 1.5; % delay b/w cs1 and cs2
consumdelay = 3;
meanITI = 60;
inhibitionstrength = -0.6; % level of dopamine response during inhibition

% ANCCR model parameter
samplingperiod = 0.2;
w = 0.5;
k = 1;
Tratio = 1.2;
alpha = 0.02;
alpha_r = 0.2;
threshold = 0.6;
minimumrate = 10^(-3);
beta = [0,0,1]';
Rtrue = [0,0,1]';
maximumjitter = 0.1;

nExp = size(meanITI,1); % number of experiments

%%
nIter = 20;
darsp = cell(length(inhibitionstrength)+1,4);
for iIter = 1:nIter
    %generate eventlog
    eventlog = simulateEventChain(numcue, 2, 4, meanITI, ...
        meanITI*3, cuecuedelay, cuerewdelay, 1, consumdelay);

    %testing with anccr model with and without DA inhibition
    for iinh = 1:length(inhibitionstrength)+1
        printf('%d iteration: %d inh condition\n',iIter,iinh);
        if iinh<=length(inhibitionstrength) % w/ inhibition 
            optolog = simulateInhibitionBlock(eventlog, 2, nan,...
                inhibitionstrength(iinh), 1, 1, numcue);
        else
            optolog = nan; % w/o inhibition 
        end
        
        [DA,ANCCR,PRC,SRC,NC,Rs] =...
            calculateANCCR(eventlog, (meanITI+cuerewdelay+consumdelay)*Tratio, alpha, k,...
            samplingperiod,w,threshold,minimumrate,beta,Rtrue,maximumjitter,optolog);

        % pool dopamine response for each event
        for ie = 1:3
           inevent = eventlog(:,1)==ie; % 1:cs1, 2:cs2, 3:reward
           darsp{iinh,ie}(:,iIter) = DA(inevent);
        end
    end
end
%% save data
dir = 'D:\OneDrive - University of California, San Francisco\figures\manuscript\dopamine_contingency\revision';

%%
close all
clr_light = [0.6,0.6,0.6; 0.6,0.6,1];
cuelist = {'CS1,','CS2'};
clr = [0 0 0; 0 0 1];

%% Fig6C middle (ANCCR)
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 4 3]);
axes('Position',axpt(5,5,2:5,1:4))
hold on;
for i= 1:2
% m = mean(darsp{2,i},2)';
% s = std(darsp{2,i},[],2)'/sqrt(nIter);
plot(1:numcue,darsp{2,i},'Color',clr_light(i,:),'LineWidth',0.35)
plot([10 30],[2 2]-i*0.2,'Color',clr(i,:),'LineWidth',1);
text(40,2-i*0.2,cuelist{i})
end
for i = 1:2
    m = mean(darsp{2,i},2)';
    plot(1:numcue,m,'Color',clr(i,:),'LineWidth',1)
end
xlabel('Trial')
ylabel('Predicted DA response')
set(gca,'FontSize',8,'Box','off','TickDir','out','LineWidth',0.35,'XTick',0:100:400,...
    'YTick',-1:2,'YLim',[-1.5 2.2],'XTickLabelRotation',45)
