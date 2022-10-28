clearvars; clc; close all;
rng(2);

%% parameter set up
% task parameters
meanIRI = 12;
numrewards = 1000;

meanITI = 3;
cuerewdelay = [3 9]; % cue-duration change, background reward, extinction
outcomedelay = [3 3];
numcue = 2000;
postrewdelay = 3;

% anccr model parameters
samplingperiod = 0.2;   
alpha_anccr = 0.02;        
alpha_r = 0.2;
w = 0.5;               
k = 1;                 
minimumrate = 10^(-3);
maximumjitter = 0.1;
beta = [0,0,1,0];
threshold = 0.6;
Tratio = 1.2;

% rpe model parameters - csc/microstimulus
alpha_rpe = 0.05;
gamma = 0.95;
lambda = 0;
statesize = 0.2;

nIter = 100;
%%
cuersp = nan(nIter,round(numcue*4/3));
for iIter = 1:nIter
    iIter
    eventlog_rr = simulateBackgroundRewards(numrewards,....
        meanIRI,3,1,0);
    [eventlog_c1,IRI_c1] = simulateEvents(repmat(numcue,1,2), [1,2], [3,nan], ...
        [1,nan], [4,nan], meanITI, meanITI*3, cuerewdelay(1), [1,0], postrewdelay);
    [eventlog_c2,IRI_c2] = simulateEvents(repmat(numcue,1,2), [1,2], [3,nan], ...
        [1,nan], [4,nan], meanITI, meanITI*3, cuerewdelay(2), [1,0], postrewdelay);
    
    eventlog = joinEventlogs(eventlog_rr,eventlog_c1,eventlog_c2);
    T = [ones(numrewards,1)*meanIRI; ones(size(eventlog_c1,1),1)*IRI_c1;...
        ones(size(eventlog_c2,1),1)*IRI_c2];
    
    [DA,ANCCR,PRC,SRC,NC] = calculateANCCR(eventlog, T*Tratio, alpha_anccr, k,...
        samplingperiod,w,threshold,minimumrate,beta,alpha_r,maximumjitter,nan,[4,3],0);
    incue = eventlog(:,1)==1;
    
    cuersp_temp = DA(incue);
    
    cuersp(iIter,:) = cuersp_temp(end-round(numcue*4/3)+1:end);
end

%%
cd('D:\OneDrive - University of California, San Francisco\figures\manuscript\dopamine_contingency\revision\data');
save('cuerewlearning_shortiri.mat','cuersp');

