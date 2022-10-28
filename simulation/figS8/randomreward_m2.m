clearvars; clc;
rng(7)

%% parameter set up
% task parameters
meanITI = 12;
numrewards = 500;

% anccr model parameters
samplingperiod = 0.2;
alpha_anccr.exponent = 0.1;
alpha_anccr.init = 0.25;
alpha_anccr.min = 0.02;        
alpha_r = 0.2;
w = 0.5;               
k = 1;                 
minimumrate = 10^(-3);
maximumjitter = 0.1;
beta = [1];
threshold = 0.6;
T = meanITI*10;

nIter = 100;
rwrsp = nan(numrewards,nIter);

%% simulation
for iiter = 1:nIter
    iiter
    eventlog = simulateBackgroundRewards(numrewards(1),meanITI,1,1,0);
    eventlog(:,2) = eventlog(:,2)+1;
    
    DA = calculateANCCR(eventlog(1:numrewards,:),T,alpha_anccr,k,samplingperiod,w,threshold,...
        minimumrate,beta,alpha_r,maximumjitter);
    eventtimeline = eventlog(1:numrewards,1);
    
    rwrsp(:,iiter) = DA(eventtimeline(:,1)==1);
end

%% figure

fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 4 3.3]);
hold on;
plot(cumsum(rwrsp,1)./repmat(sum(rwrsp,1),numrewards,1),'Color',[0.6 0.6 0.6],'LineWidth',0.35);
plot(mean(cumsum(rwrsp,1)./repmat(sum(rwrsp,1),numrewards,1),2),'k','LineWidth',1);
ylim([-1.3 1.5])
xlim([0 numrewards])
xlabel('Trial');
ylabel({'Normalized'; 'cumsum (predicted DA)'})
set(gca,'Box','off','TickDir','out','FontSize',8,'LineWidth',0.35,'XTick',...
    [0 200 400],'YTick',-1:1:1);
cd('D:\OneDrive - University of California, San Francisco\figures\manuscript\dopamine_contingency\revision');
print(fHandle,'-depsc','-painters','randomreward_m2.ai')

%%
% fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 4 3.3]);
% scatter(1:200,rwrsp(1:200,1),3,'k','filled')
% hold on;
% beta = glmfit(1:200,rwrsp(1:200,1));
% plot([1 200],[1 200]*beta(2)+beta(1),'r');
% [r,p]= corr([1:200]',rwrsp(1:200,1));
% xlabel('Trial');
% ylabel('Predicted DA response');
% text(120,-0.4,['r = ',num2str(round(r*1000)/1000)],'FontSize',5);
% text(120,-0.6,['p = ',num2str(p)],'FontSize',5);
% set(gca,'Box','off','TickDir','out','FontSize',5,'LineWidth',0.35,'XTick',0:100:200,'YTick',-1:2);
% print(fHandle,'-dtiff','-r600','randomreward_anccr.tif')