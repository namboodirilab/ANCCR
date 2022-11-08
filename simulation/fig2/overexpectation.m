clearvars; clc; 

rng(2)

% Task parameter
meanITI = [50, 50];
maxITI = meanITI*3;
cuerewdelay = [0.5, 0.5]; 
num_training = [500, 500];
num_paired = 500;
num_test = 3;

% Model parameter
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

T = Tratio*(meanITI(1)+cuerewdelay(1));

nIter = 20;
% Mark new idxs for giving CS1 only w/o CS2 or reward
cs1_testidxs = [num_training(1)+[1:3], num_training(1)+num_paired+ [4:6]];

for iIter = 1:nIter
    % Generate 200 instances (each) of C1 ->R, C2 -> R
    eventlog_training = simulateEvents(num_training, [1, 2], [3, 3], [1, 1], nan,...
        meanITI, maxITI, cuerewdelay, [1, 1], [0, 0]);
    % Generate 3 test C1 instances
    eventlog_test1 = simulateEvents(num_test, 1, 3, 1, nan,...
        meanITI, maxITI, cuerewdelay, 0, 0);
    eventlog_test1(:,2) = eventlog_test1(:,2)+eventlog_training(end,2);
    % Generate 200 instances of C1 + C2 -> R
    eventlog_paired = simulateEventChain(num_paired, 2, nan, meanITI(1), ...
        maxITI(1), 0, cuerewdelay(1), 1, 0);
    eventlog_paired(:,2) = eventlog_paired(:,2)+eventlog_test1(end,2);
    % Generate 3 more test C1 instances
    eventlog_test2 = simulateEvents(num_test, 1, 3, 1, nan,...
        meanITI, maxITI, cuerewdelay, 0, 0);
    eventlog_test2(:,2) = eventlog_test2(:,2)+eventlog_paired(end,2);
    % Concat. all generated events
    eventlog = [eventlog_training; eventlog_test1; eventlog_paired; 
        eventlog_test2];
    % Calculate model values
    [~,~,~,SRC,~,R] =...
        calculateANCCR(eventlog,T,alpha,k,samplingperiod,w,threshold,minimumrate,...
        beta,alpha_r,maximumjitter,nan);
    for ie = 1:3
        % Save SRC, Rs across iterations for value calculation
        SRCs{ie}(iIter,:) = squeeze(SRC(ie,3,eventlog(:,1)==ie));
        Rs{ie}(iIter,:) = squeeze(R(ie,3,eventlog(:,1)==ie));
    end
end

% Calculate value
Beh = cellfun(@(x,y) x.*y,SRCs,Rs,'UniformOutput',false);

%% Plotting SRC*R (rough value of cue) response
clr = {'k','r',;[0.6 0.6 0.6],[1 0.6 0.6]};

beta= 5;
temperature = 1/beta;
cost = -0.3;

% Calculate action probabilities
prob = cell(2,1);
for itest = 1:length(cs1_testidxs)
    q = [zeros(size(Beh{1},1),1), Beh{1}(:,cs1_testidxs(itest))+cost];
    prob{1}(:,itest) = exp(q(:,2)/temperature)./sum(exp(q/temperature),2);
end

% Plot action probabilities
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 2.5 3.5]);
axes('Position',axpt(5,5,2:5,1:4)) 
hold on;
data_training = mean(prob{1}(:,1:3),2);
data_paired = mean(prob{1}(:,4:end),2);
bar(1,mean(data_training),0.5,'FaceColor',clr{2,1},'EdgeColor','none');
bar(2,mean(data_paired),0.5,'FaceColor',clr{2,2},'EdgeColor','none');
scatter(rand(nIter,1)*0.4-0.2+1,data_training,5,clr{1,1},'filled')
scatter(rand(nIter,1)*0.4-0.2+2,data_paired,5,clr{1,2},'filled')
errorbar(1,mean(data_training),std(data_training)/sqrt(nIter),clr{1,1})
errorbar(2,mean(data_paired),std(data_paired)/sqrt(nIter),clr{1,2})

set(gca,'Box','off','TickDir','out','FontSize',8,'LineWidth',0.35,...
    'XTick',[1,2],'XTickLabel',{'Prior to compound training', ...
    'After compound training'},'XTickLAbelRotation',45);
xlim([0.5 2.5])
ylim([0.5 1])
ylabel('P(lick|CS1)')

savefig(fHandle,'overexpectation.fig')
print(fHandle,'-depsc','-painters','overexpectation.ai')
save('overexpectation.mat', 'prob', 'Beh')

[significant_binary, p_value, ~, stats] = ttest(data_training, data_paired);

close all