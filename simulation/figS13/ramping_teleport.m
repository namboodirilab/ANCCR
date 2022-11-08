% simulate teleport task used in Kim et al., 2020 with ANCCR

clearvars; close all; clc;
rng(2)

% task parameter
numcue = 2000;
cuerewdelay = 1; % delay from cs2 to reward
cuecuedelay = 1; % delay b/w cs1 and cs2
consumdelay = 3;
meanITI = 6;
speed = [0.5,1,2];

% model parameter
samplingperiod = 0.2;
w = 0.5;
k = 1;
Tratio = 1.2;
alpha = 0.02;
alpha_r = 0.2;
threshold = 0.4;
minimumrate = 10^(-3);
beta = [zeros(1,8),1]';
maximumjitter = 0.1;

nExp = size(meanITI,1); % number of experiments

%%

darsp = cell(4,1);
nIter = 100;
for iiter = 1:nIter
    iiter
    %% generate eventlog
    % training before incorporating teleport trials
    eventlog_pre = simulateEventChain(numcue, 8, nan, meanITI, ...
        meanITI*3, cuecuedelay, cuerewdelay, 1, consumdelay);
    % test trials including teleport trials
    eventlog_post = simulateEventChain(numcue, 8, nan, meanITI, ...
        meanITI*3, cuecuedelay, cuerewdelay, 1, consumdelay);

    % in 15% of total trials, animal was teleported either 
    % from 1s to 2s, 3s to 4s, or 5s to 6s from first cue onset
    testidx = randsample(find(eventlog_post(:,1)==1),round(numcue*0.15));
    testidx = cellfun(@(x) sort(testidx(round(numcue*0.03*x/3)+1:round(numcue*0.03*(x+1)/3))),...
        num2cell(0:2),'UniformOutput',false);
    [~,testtrial] = cellfun(@(x) ismember(x,find(eventlog_post(:,1)==1)),testidx,'UniformOutput',false);
    testtrial = cellfun(@(x) x+numcue,testtrial,'UniformOutput',false);
   
    for itest = 1:3
        eventlog_post(testidx{itest}+1+(itest-1)*2,:) = nan;
        eventlog_post([itest*2:8]+testidx{itest},2) =...
            eventlog_post([itest*2:8]+testidx{itest},2)-cuecuedelay;
    end
    eventlog_post = rmmissing(eventlog_post);

    % control trial w/o teleport
    ctrltrial = numcue+1:numcue*2;
    ctrltrial = ctrltrial(~ismember(ctrltrial,cell2mat(testtrial)));
    
    eventlog = joinEventlogs(eventlog_pre,eventlog_post);

    %% simulate
    [DA,ANCCR,PRC,SRC,NC] = calculateANCCR(eventlog, (meanITI+cuerewdelay+consumdelay)*Tratio, alpha, k,...
        samplingperiod,w,threshold,minimumrate,beta,alpha_r,maximumjitter);

    % pool dopamine response
    incue = find(eventlog(:,1)==1);
    for itest = 1:3
        intest = incue(testtrial{itest})+1+(itest-1)*2;
        darsp{itest}(iiter,:) = mean(DA(incue(testtrial{itest})+[0:6]),1);
    end
    darsp{4}(iiter,:) = mean(DA(incue(ctrltrial)+[0:7]),1);
end

%% save data
dir = 'D:\OneDrive - University of California, San Francisco\figures\manuscript\dopamine_contingency\revision\';
cd(dir)
save('data\ramping_teleport.mat','darsp','threshold');

%% FigS13B
ct = cbrewer('seq','YlOrRd',5);
ct = [flip(ct([2,4:5],:));0 0 0];

close all
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 3.2 3.2]);
hold on;
for i = 1:3
   bar(i,nanmean(darsp{i}(:,2*i)./darsp{4}(:,2*i)),0.8,'FaceColor',ct(i,:),'EdgeColor','none');
   errorbar(i,nanmean(darsp{i}(:,2*i)./darsp{4}(:,2*i)),...
       nanstd(darsp{i}(:,2*i)./darsp{4}(:,2*i))/sqrt(size(darsp{i},1)),'k','CapSize',3);
   [~,p,~,stat] = ttest(darsp{i}(:,2*i),darsp{4}(:,2*i))
end
bar(4,1,0.8,'FaceColor',ct(4,:),'EdgeColor','none');
plot([0 5],[1 1],'k--');
xlim([0.25 4.75]);
ylim([0.5 1.5]);
ylabel({'Normalized';'Predicted DA'});
set(gca,'Box','off','TickDir','out','FontSize',8,'LineWidth',0.35,...
    'XTick',1:4,'XTickLabel',{'T1';'T2';'T3';'Stand.'},'XTickLabelRotation',45,'YTick',0.5:0.5:1.5);
