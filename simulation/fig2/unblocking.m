clearvars; clc; 

rng(5)

% task parameter
meanITI = 100;
maxITI = meanITI*3;
cuerewdelay = 0.5; 
numcs1only = 4000;
numcs1cs2 = 1000;

% model parameter
samplingperiod = 0.2;
alpha = 0.02;
alpha_r = 10*alpha;
w = 0.5;
k = 0.01;
Tratio = 1.2;
minimumrate = 10^(-3);
threshold = 0.6;
maximumjitter = 0.1;
beta = [0, 0, 1]';
Rtrue = [0, 0, 1]';

exact_mean_or_not = 0;

nIter = 20;
testtrial = numcs1only+[300:302]; % giving CS2 only w/o CS1
%%
[SRCs,Rs] = deal(cell(3,2));
for iIter = 1:nIter
    iIter
    %% w/ blocking: C1-R --> C1&C2-R
    % pre-blocking
    [eventlog,IRI] = simulateEvents(numcs1only, 1, 3, 1, nan,...
        meanITI, meanITI*3, cuerewdelay, 1, 0);
    
    eventlog_temp = simulateEventChain(numcs1cs2, 2, nan, meanITI, ...
        meanITI*3, 0, cuerewdelay, 1, 0);
    eventlog_temp(:,2) = eventlog_temp(:,2)+eventlog(end,2);
    eventlog = [eventlog;eventlog_temp];

    % probe test
    cs1idx = find(eventlog(:,1)==1);
    rwidx = find(eventlog(:,1)==3);
    eventlog([cs1idx(testtrial),rwidx(testtrial)],:) = [];

    for iAct = 1:2
        if iAct==1
            optolog = nan;
        else
            optolog = simulateInhibitionBlock(eventlog, 3, nan, 1, 1, numcs1only, nan);
        end
        
        [DA,~,PRC,SRC,NC,R] =...
            calculateANCCR(eventlog,IRI*Tratio,alpha,k,samplingperiod,w,threshold,minimumrate,...
            beta,alpha_r,maximumjitter,optolog);
        for ie = 1:3
            SRCs{ie,iAct}(iIter,:) = squeeze(SRC(ie,3,eventlog(:,1)==ie));
            Rs{ie,iAct}(iIter,:) = squeeze(R(ie,3,eventlog(:,1)==ie));
        end
    end
end

%%
Beh = cellfun(@(x,y) x.*y,SRCs,Rs,'UniformOutput',false);

%% plotting behavior 
%   P(a,s) = softmax(Q(a,s))*Indicator(s)
%   possible (a,s): no lick, lick after cs1, lick after cs2 
%   Q(nolick,cs) = 0;
%   Q(lick,cs) = SRC(cs,reward)*R(cs)
%   Indicator(s) = whether or not s is when meaningful causal target is given
%%
beta= 5;
temperature = 1/beta;
cost = -0.3;

prob = cell(2,1);
for i = 1:2
for itest = 1:length(testtrial)
q = [zeros(size(Beh{2,i},1),1), Beh{2,i}(:,testtrial(itest)-numcs1only)+cost];
prob{i}(:,itest) = exp(q(:,2)/temperature)./sum(exp(q/temperature),2);
end
end

%%
clr = {[0 0 0],[0 0 1];[0.6 0.6 0.6],[0.6 0.6 1]};
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 2.5 3.5]);
axes('Position',axpt(5,5,2:5,1:4)) 
hold on;
for i = 1:2
data = mean(prob{i}(:,1:3),2);
bar(i,mean(data),0.5,'FaceColor',clr{2,i},'EdgeColor','k');
scatter(rand(nIter,1)*0.4-0.2+i,data,5,clr{1,i},'filled')
errorbar(i,mean(data),std(data)/sqrt(nIter),'Color','k')
end
set(gca,'XTick',[1,2],'XTickLabel',{'w/o act','w/ act'},...
    'XTickLAbelRotation',45,'YTick',0:0.5:1,'TickDir','out','LineWidth',0.35);
xlim([0.5 2.5])
ylim([0 1.1])
ylabel('P(lick|CS2)')

dir = 'D:\OneDrive - University of California, San Francisco\dopamine contingency\figures\figures for r3';
cd(dir);
print(fHandle,'-depsc','-painters','unblocking.ai')


%%

[~,p,~,stat] = ttest2(mean(prob{1}(:,1:3),2),mean(prob{2}(:,1:3),2));
