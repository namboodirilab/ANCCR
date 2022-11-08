% simulate pavlovian conditioning task w/o pre-conditioning reward exposure
% (simulate similar condition used in Coddington and Dudman, 2018)

clearvars; clc; close all;

rng(2)

% task parameter
meanITI = 28;
cuerewdelay = 1.5;
outcomedelay = 0;
numcue = 2000;

% model parameter
samplingperiod = 0.2;
alpha.exponent = 0.1;
alpha.init = 0.25;
alpha.min = 0.02;  
alpha_r = 0.2;
w = 0.5;
k = 0.05;
minimumrate = 10^(-3);
maximumjitter = 0.1;
beta = [0,1];
threshold = 0.6;
Tratio = 1.2;

nIter = 100;
rsprw = nan(numcue,nIter);
%% experiment set up

for iiter = 1:nIter
    iiter
    % generate eventlog
    [eventlog,IRI] = simulateEvents(numcue, 1, 2, ...
        1, nan, meanITI, meanITI*3, cuerewdelay, 1,outcomedelay);

    %simulate
    DA = calculateANCCR(eventlog, IRI*Tratio, alpha, k,...
        samplingperiod,w,threshold,minimumrate,beta,alpha_r,maximumjitter);

    %pool reward response
    rsprw(:,iiter) = DA(eventlog(:,1)==2);

end

%% save data
dir = 'D:\OneDrive - University of California, San Francisco\figures\manuscript\dopamine_contingency\revision';

%% FigS13E
m = mean(rsprw,2)';
s = std(rsprw,[],2)'/sqrt(nIter);
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 4 3.2]);
hold on;
fill([1:numcue flip(1:numcue)],[m+s flip(m-s)],[0.6 0.6 0.6],'EdgeColor','none');
plot(1:2000,m,'k');
plot([0 500],[0 0],'k:','LineWidth',0.35);
set(gca,'Box','off','TickDir','out','FontSize',8,'XTick',0:250:500,'YTick',0:0.5:1);
ylim([-0.1 1]);
xlim([0 500]);
xlabel('Trial')
ylabel({'Predicted'; 'reward response'});

