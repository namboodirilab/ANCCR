% simulate sequential learning task with gradually increasing inhibition
% level using microstimulus model

clearvars; close all; clc;
rng(2)

% task parameter
numcue = 800;
cuerewdelay = 3;
cuecuedelay = 1.5;
consumdelay = 3;
meanITI = 60;
inhibitionstrength = -0.6;

% model parameter
alpha = 0.05; % learning rate
gamma = 0.95; % discount factor
lambda = 0.95;
numstimulus = 20; % number of microstimulus
d = 0.99; % memory decay constant
statesize = 0.2;
sigma = 0.08; % width of basis function

%%

nIter = 20;
darsp = cell(2,3);
for iIter = 1:nIter
    % generate eventlog
    eventlog = simulateEventChain(numcue, 2, 4, meanITI, ...
        meanITI*3, cuecuedelay, cuerewdelay, 1, consumdelay);

    % inhibition magnitude starts from 0 and linearly increases to -0.6
    inhibitionlog = [eventlog(eventlog(:,1)==2,2),...
        eventlog(eventlog(:,1)==2,2)+cuecuedelay,...
        zeros(numcue,1), ones(numcue,1)*inhibitionstrength]; 

    % simulate 
    [rpetimeline,valuetimeline,eventtimeline,inhibitiontimeline] =...
        simulatemicrostimulus(eventlog,3,numstimulus,statesize,alpha,gamma,lambda,sigma,d,inhibitionlog);
    
    for ie = 1:3 % 1:cs1, 2:cs2, 3:reward
        inevent = eventtimeline(:,1)==ie;
        darsp{iinh,ie}(:,iIter) = rpetimeline(inevent);
    end
end

%% FigS12D
close all
clr ={'r','k'};
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 2 3]);
hold on;
iinh = 2;
    data = cell2mat(cellfun(@(x) mean(x(end-100:end,1:nIter),1),darsp(iinh,1),'UniformOutput',false))';
        scatter(rand(nIter,1),data,2,[0.6 0.6 0.6],'filled');
        errorbar(0.5,mean(data),std(data)/sqrt(nIter),'LineWidth',0.5,'Color','k');
plot([-0.5 1.5],[0 0],'k:','LineWidth',0.35);
ylim([-4 1])
xlim([-0.5 1.5]);
ylabel({'Predicted';'CS1 response'});
title({'Model 1';'microstimluus'})
set(gca,'Box','off','TickDir','out','FontSize',8,'LineWidth',0.35,'XTick',0.5,...
    'XTickLabel',[],'YTick',-4:2:0);

