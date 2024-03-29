% simulate speed change task used in Kim et al., 2020 with ANCCR

clearvars; close all; clc;
rng(2)

% task parameter
numcue = 2000;
cuerewdelay = 1; % delay from cs2 to reward
cuecuedelay = 1; % delay b/w cs1 and cs2
consumdelay = 3;
meanITI = 6;
speed = [0.5,2]; % scale of speed: 0.5,faster; 2,slower
IRI = cuerewdelay*8+meanITI+consumdelay;

% model parameter
samplingperiod = 0.2;
w = 0.5;
k = 1;
Tratio = [1.2;0.5];
exponent = 0.01;
alpha = 0.02;
alpha_r = 0.2;
threshold = 0.5;
minimumrate = 10^(-3);
beta = [zeros(1,8),1]';
maximumjitter = 0.1;

nExp = size(meanITI,1); % number of experiments

%%
darsp = cell(3,1);
nIter = 100;
for iiter = 1:nIter
    iiter  
    %% generate eventlog
    % training before incorporating speed change trials
    eventlog_pre = simulateEventChain(numcue, 8, nan, meanITI, ...
        meanITI*3, cuecuedelay, cuerewdelay, 1, consumdelay);
    % test trials including speed change trials
    eventlog_post = simulateEventChain(numcue, 8, nan, meanITI, ...
        meanITI*3, cuecuedelay, cuerewdelay, 1, consumdelay);

    % in 20% of total trials, speed was manipulated to X2 or X0.5
    testidx = randsample(find(eventlog_post(:,1)==1),round(numcue*0.2));
    testidx = cellfun(@(x) sort(testidx(round(numcue*0.2*x/2)+1:round(numcue*0.2*(x+1)/2))),...
        num2cell(0:1),'UniformOutput',false);
    [~,testtrial] = cellfun(@(x) ismember(x,find(eventlog_post(:,1)==1)),testidx,'UniformOutput',false);
    testtrial = cellfun(@(x) x+numcue,testtrial,'UniformOutput',false);

    for itest = 1:2
        for i = 1:length(testtrial{itest})
            deltat = 8*(speed(itest)-1)/cuecuedelay;
            eventlog_post([0:8]+testidx{itest}(i),2) =...
                eventlog_post(testidx{itest}(i),2)+[0:8]*speed(itest)/cuecuedelay;
            eventlog_post(9+testidx{itest}(i):end,2) =...
                eventlog_post(9+testidx{itest}(i):end,2)+deltat;
        end
    end
    
    % control trial w/o speed change 
    ctrltrial = numcue+1:numcue*2;
    ctrltrial = ctrltrial(~ismember(ctrltrial,cell2mat(testtrial)));
   
    %% simulate 
    eventlog = joinEventlogs(eventlog_pre,eventlog_post);
    
    cueidx = find(eventlog(:,1)==1);
    speedlog = ones(size(eventlog,1),1);
    speedlog(cueidx(testtrial{1})+[0:8]) = 1/speed(1);
    speedlog(cueidx(testtrial{2})+[0:8]) = 1/speed(2);
    
    % assume change in T
    T = [repmat(IRI*Tratio(1),size(eventlog_pre,1),1);...
        exp(-exponent*[1:size(eventlog_post,1)]')*(IRI*Tratio(1)-IRI*Tratio(2))+IRI*Tratio(2)];
    [DA,ANCCR,PRC,SRC,NC] = calculateANCCR(eventlog, T, alpha, k,...
        samplingperiod,w,threshold,minimumrate,beta,alpha_r,maximumjitter); % now speed is multiplied directly to NC
    
    % pool dopamine response
    incue = find(eventlog(:,1)==1);
    for itest = 1:2
        intest = incue(testtrial{itest})+1+(itest-1)*2;
        darsp{itest}(iiter,:) = mean(DA(incue(testtrial{itest})+[0:7]),1);
        darsp{itest}(iiter,2:end) = darsp{itest}(iiter,2:end)*(1/speed(itest));
    end
    darsp{3}(iiter,:) = mean(DA(incue(ctrltrial)+[0:7]),1); % 1:slow, 2:fast, 3:standard
end

%% FigS13D
dir = 'D:\OneDrive - UCSF\dopamine contingency\erratum\';
close all
ct = cbrewer('seq','YlOrRd',5);
ct = [ct([2,5],:);0 0 0];
ratio = [speed,1];

x = [3,1,2];
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 2.6 3.2]);
hold on;
for i = 1:3
    data(:,i) = darsp{i}(:,end)./darsp{3}(:,end);
    bar(x(i),nanmean(data(:,i)),0.8,'FaceColor',ct(i,:),'EdgeColor','none');
    errorbar(x(i),nanmean(data(:,i)),nanstd(data(:,i))/sqrt(size(data,1)),'k','CapSize',3);
    [~,p,~,stat] = ttest(darsp{i}(:,end),darsp{3}(:,end));
end
plot([0 4],[1 1],'k--','LineWidth',0.5);
xlim([0.25 3.75]);
ylim([0 2]);
ylabel({'Normalized';'predicted DA'});
set(gca,'Box','off','TickDir','out','FontSize',8,'LineWidth',0.35,...
    'XTick',1:4,'XTickLabel',{'Slow';'Stand.';'Fast'},'XTickLabelRotation',45,'YTick',0:1:2);
print(fHandle,'-depsc',[dir,'ramping_speed.ai']);