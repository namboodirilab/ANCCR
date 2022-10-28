clearvars; clc; 

rng(2)

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
testtrial = [numcs1only+[300:302],numcs1cs2+numcs1only-2:numcs1cs2+numcs1only]; % giving CS2 only w/o CS1
%%
for iIter = 1:nIter
    %% w/ blocking: C1-R --> C1&C2-R
    % pre-blocking
    eventlog = simulateEvents(numcs1only, 1, 3, 1, nan,...
        meanITI, meanITI*3, cuerewdelay, 1, 0);
    
    eventlog_temp = simulateEventChain(numcs1cs2, 2, nan, meanITI, ...
        meanITI*3, 0, cuerewdelay, 1, 0);
    eventlog_temp(:,2) = eventlog_temp(:,2)+eventlog(end,2);
    eventlog = [eventlog;eventlog_temp];

    % probe test
    cs1idx = find(eventlog(:,1)==1);
    rwidx = find(eventlog(:,1)==3);
    eventlog([cs1idx(testtrial),rwidx(testtrial)],:) = [];
    T = Tratio*(meanITI+cuerewdelay);

    [DA,~,PRC,SRC,NC,R] =...
        calculateANCCR(eventlog,T,alpha,k,samplingperiod,w,threshold,minimumrate,...
        beta,alpha_r,maximumjitter,nan);
    for ie = 1:3
        SRCs{ie,1}(iIter,:) = squeeze(SRC(ie,3,eventlog(:,1)==ie));
        Rs{ie,1}(iIter,:) = squeeze(R(ie,3,eventlog(:,1)==ie));
    end

    %% w/o blocking: C1-R 
    eventlog = simulateEvents(numcs1cs2, 1, 2, 1, nan,...
        meanITI, meanITI*3, cuerewdelay, 1, 0);
    
    [DA,~,PRC,SRC,~,R] =...
        calculateANCCR(eventlog,T,alpha,k,samplingperiod,w,threshold,minimumrate,...
        beta(2:end),alpha_r,maximumjitter,nan);
    for ie = 2:3
        SRCs{ie,2}(iIter,:) = squeeze(SRC(ie-1,3-1,eventlog(:,1)==ie-1));
        Rs{ie,2}(iIter,:) = squeeze(R(ie-1,3-1,eventlog(:,1)==ie-1));
    end
end

%%
Beh = cellfun(@(x,y) x.*y,SRCs,Rs,'UniformOutput',false);

%% plotting SRC*R (rough value of cue) response
clr = {'k','r',;[0.6 0.6 0.6],[1 0.6 0.6]};
close all
figure
subplot(2,1,1)
hold on;
for i = 3:-1:1
    if i==2
        x = [1:numcs1cs2]+numcs1only;
    else
        x = 1:numcs1only+numcs1cs2;
        x(testtrial) = [];
    end
plot(x,mean(Beh{i,1},1));
end
% x = [repmat(testtrial,2,1);NaN,NaN];
% y = repmat([-1;1.5;NaN],1,2);
% plot(x(:),y(:),'k');
plot([numcs1only numcs1only],[-1 1.5],'k:')
legend({'R';'C2';'C1'},'Location','southwest')
ylabel('SRC(x,R)*R(x)')
ylim([-0.5 1.5])

subplot(2,1,2)
hold on;
plot([1:numcs1cs2],mean(Beh{2,1},1),clr{1,1});
plot([1:numcs1cs2],mean(Beh{2,2},1),clr{1,2});
legend({'blocking';'non-blocking'},'Location','northwest')
ylabel('SRC(C2,R)*R(C2)')

%%
figure
for i = 1:2
    hold on;
    data = [mean(Beh{2,i}(:,testtrial(1:3)-numcs1only),2),...
        mean(Beh{2,i}(:,testtrial(4:6)-numcs1only),2)];
    plot(1:2,data,'Color',clr{2,i},'LineWidth',0.35);
    plot(1:2,mean(data),...
        'Color',clr{1,i},'LineWidth',1);
end
ylabel('SRC(C2,R)*R(C2)')
set(gca,'XTick',1:2,'XTickLabel',{'early','late'},'XTickLabelRotation',45)
xlabel('Trial')
title('Probe trial')
xlim([0.5 2.5])


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

a = prob{1}; b = prob{2};
prob{1} = b; prob{2} = a;

fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 2.5 3.5]);
axes('Position',axpt(5,5,2:5,1:4)) 
hold on;
for i = 1:2
    data = mean(prob{i}(:,1:3),2);
    bar(i,mean(data),0.5,'FaceColor',clr{2,i},'EdgeColor','none');
    scatter(rand(nIter,1)*0.4-0.2+i,data,5,clr{1,i},'filled')
    errorbar(i,mean(data),std(data)/sqrt(nIter),clr{1,i})
end
set(gca,'Box','off','TickDir','out','FontSize',8,'LineWidth',0.35,...
    'XTick',[1,2],'XTickLabel',{'Control','Blocking'},'XTickLAbelRotation',45);
xlim([0.5 2.5])
ylabel('P(lick|CS1)')

savefig(fHandle,'blocking.fig')
print(fHandle,'-depsc','-painters','blocking.ai')
save('blocking.mat', 'prob', 'Beh')

control = mean(prob{1}(:,1:3),2);
block = mean(prob{2}(:,1:3),2);
[significant_binary, p_value, ~, stats] = ttest(control, block);

close all

%%
[~,p,~,stat] = ttest2(mean(prob{1}(:,1:3),2),mean(prob{2}(:,1:3),2))
