%DA signals +/- prediction errors to US and CS, 
% in proportion to reward magnitude and/or probability of reward.
clearvars; clc; 

rng(2)
%% Task setup
% task parameter
cuerewdelay = [1, 1, 1];
postrewdelay = [1, 1, 1] ;
numcue = [50, 900, 50];
rew_mags = [1; 5; 10];
rew_probs = [1; 1; 1];
meanITI = [30 30 30];
maxITI = meanITI*3;
IRI = meanITI + cuerewdelay;

% model parameter
samplingperiod = 0.2;
alpha = 0.02;
alpha_r = alpha*10;
w = 0.5;
k = 1;
Tratio = 1.2;
minimumrate = 10^(-3);
threshold = 0.6;
maximumjitter = 0.1;
beta = [0; 1];
Rtrue = [0; 1];
nIter = 20;
da_means = nan(3, nIter);
da_stds = nan(3, nIter);

%% Run simulation
%Simulate for cue reward pairs with three different magnitudes
for iIter = 1:nIter
    [eventlog] = simulateEvents(numcue, [1,1,1], [2,2,2], rew_mags, nan,...
    meanITI, maxITI, cuerewdelay, rew_probs, postrewdelay);

    [DA,ANCCR,PRC,SRC,NC,Rs] = calculateANCCR(eventlog, IRI(1)*Tratio, alpha, k, ...
    samplingperiod,w,threshold,minimumrate,beta,alpha_r, ...
    maximumjitter,nan,nan);

    tens = find(eventlog(:,3) == 10); fives = find(eventlog(:,3) == 5);
    ones = find(eventlog(:,3) == 1);
    da_10 = DA(tens); da_5 = DA(fives); da_1 = DA(ones);
    da_means(:, iIter) = [mean(da_1), mean(da_5), mean(da_10)];
end

%%
close all
dir = 'D:\OneDrive - University of California, San Francisco\figures\manuscript\dopamine_contingency\revision';

fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 2.5 3.5]);
axes('Position',axpt(5,5,2:5,1:4)) 
hold on;
plot(1:3,da_means,'Color',[0.6 0.6 0.6],'LineWidth',0.35);
errorbar(1:3,mean(da_means,2),std(da_means,[],2)/sqrt(nIter),'k','LineWidth',0.5);
plot([0 4],[0 0],'k:','LineWidth',0.35);
set(gca,'Box','off','TickDir','out','FontSize',8,'LineWidth',0.35,...
    'XTick',1:3,'XTickLabel',rew_mags,...
    'YTick',-4:4:8,'YLim',[-4 8],'XLim',[0.5 3.5])
xlabel('Reward magnitude')
ylabel('Predicted DA response')
print(fHandle,'-depsc','-painters',[dir,'\rew_mag.ai']);

% 
% f = figure; 
% hold on;
% for i = 1:nIter
%     errorbar([1, 2, 3], da_means(:, i), da_stds(:, i), ':k',"LineWidth",0.1)
% end
% 
% errorbar([1, 2, 3], mean(da_means,2), mean(da_stds,2),'k',"LineWidth", 2)
% figure; hold on;
% for i = 1:nIter
%     errorbar([1, 2, 3], da_means(:, i), da_stds(:, i), ':k',"LineWidth",0.5)
% end
% 
% errorbar([1, 2, 3], mean(da_means,2), mean(da_stds,2),'k',"LineWidth",1.5)
% xticks([1, 2, 3])
% xticklabels({'r=10', 'r=5', 'r=1'})
% axis padded
% ylabel('DA(r)')
% title('Dopamine response to trial-by-trial changes in reward magnitude')
% 
% savefig(f,'rewardmag_trialbytrial.fig')
% print(f,'-depsc','-painters','rewardmag_trialbytrial.ai')
% save('rewardmag_trialbytrial.mat', 'da_means', 'da_stds')
% 
% close all




