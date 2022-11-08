% Simulate blocking experiment.
clearvars; clc; 
rng(2) % Use the same seed

% Set task parameters
meanITI = 100;
maxITI = meanITI*3;
cuerewdelay = 0.5; 
numcs1only = 4000;
numcs1cs2 = 1000;

% Model parameters
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
%% Simulate
for iIter = 1:nIter
    %% w/ blocking: C1-R --> C1&C2-R
    % Pre-blocking - generate eventlog with only CS1->R
    eventlog_preblock = simulateEvents(numcs1only, 1, 3, 1, nan,...
        meanITI, meanITI*3, cuerewdelay, 1, 0);
    
    % Generate eventlog for CS1+CS2 pairings
    eventlog_paired = simulateEventChain(numcs1cs2, 2, nan, meanITI, ...
        meanITI*3, 0, cuerewdelay, 1, 0);
    % Shift eventlog for CS1+CS2 pairings by the last timestamp from 
    % eventlog_preblock
    eventlog_paired(:,2) = eventlog_paired(:,2) + eventlog_preblock(end,2);
    % Concatenate preblock and paired eventlogs
    eventlog = [eventlog_preblock; eventlog_paired];

    % Probe test
    cs1idx = find(eventlog(:,1)==1);
    rwidx = find(eventlog(:,1)==3);
    % Remove rewards for three CS1 instances for testing
    eventlog([cs1idx(testtrial),rwidx(testtrial)],:) = [];
    % Set decay time constant
    T = Tratio*(meanITI+cuerewdelay); 
    % Calculate model values
    [DA,~,PRC,SRC,NC,R] =...
        calculateANCCR(eventlog,T,alpha,k,samplingperiod,w,threshold,minimumrate,...
        beta,alpha_r,maximumjitter,nan);
    % Save predecessor representation and learned reward values for all
    % events relative to reward
    for ie = 1:3
        SRCs{ie,1}(iIter,:) = squeeze(SRC(ie,3,eventlog(:,1)==ie));
        Rs{ie,1}(iIter,:) = squeeze(R(ie,3,eventlog(:,1)==ie));
    end

    %% w/o blocking: C1-R (control)
    % Generate eventlog with only CS2->R (labels are offset by -1 relative
    % to previous block)
    eventlog = simulateEvents(numcs1cs2, 1, 2, 1, nan,...
        meanITI, meanITI*3, cuerewdelay, 1, 0);
    % Calculate model values
    [DA,~,PRC,SRC,~,R] =...
        calculateANCCR(eventlog,T,alpha,k,samplingperiod,w,threshold,minimumrate,...
        beta(2:end),alpha_r,maximumjitter,nan);
    % Save predecessor representation and learned reward values for all
    % events relative to reward
    for ie = 2:3
        SRCs{ie,2}(iIter,:) = squeeze(SRC(ie-1,3-1,eventlog(:,1)==ie-1));
        Rs{ie,2}(iIter,:) = squeeze(R(ie-1,3-1,eventlog(:,1)==ie-1));
    end
end

%% Calculate value (SRC*R)
Beh = cellfun(@(x,y) x.*y,SRCs,Rs,'UniformOutput',false);

% Plotting SRC*R (rough value of cue) response (not used in manuscript)
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


%% Use value estimate to calculate and plot action probabilities
beta= 5;
temperature = 1/beta;
cost = -0.3;

prob = cell(2,1);
for i = 1:2 % For C2 in blocking case i = 1, C2 alone i = 2
for itest = 1:length(testtrial) % For each test instance (in probe test)
% Value
q = [zeros(size(Beh{2,i},1),1), Beh{2,i}(:,testtrial(itest)-numcs1only)+cost];
% Apply softmax to find p(lick|C2)
prob{i}(:,itest) = exp(q(:,2)/temperature)./sum(exp(q/temperature),2);
end
end

% Reverse order of cell content s/t that control is plotted left
a = prob{1}; b = prob{2};
prob{1} = b; prob{2} = a;

% Plot lick probabilities for control and test case
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

close all

%% T-test
control = mean(prob{1}(:,1:3),2);
block = mean(prob{2}(:,1:3),2);
[~,p,~,stat] = ttest2(control,block);
