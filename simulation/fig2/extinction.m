% Show extinction of DA response (also behavior)
clearvars; clc; 
rng(2)

%% Initialization
% Task parameters
meanITI = 30;
maxITI = meanITI*3;
intercuedelay = 1;
cuerewdelay = 1;
postrewdelay = 1;
n_cues = 400;
inh_start = 300;

% model parameters
samplingperiod = 0.2;
alpha = 0.02;
alpha_r = alpha*10;
w = 0.5;
k = 1;
Tratio = 1.2;
T = meanITI*Tratio;
minimumrate = 10^(-3);
threshold = 0.8;
maximumjitter = 0.1;
beta = [0; 1];
exact_mean_or_not = 0;
nIter = 20;

cost = -0.5;

SRC_all_iters = nan(length(beta), length(beta), n_cues, nIter);
Rs_all_iters = nan(length(beta), length(beta), n_cues, nIter);
q_src_all_iters = nan(length(beta), n_cues, nIter);

%% Simulation
for iIter = 1:nIter
    % Simulate all cue and reward times
    [eventlog] = simulateEvents(n_cues, 1, 2, 1, nan, ...
                meanITI, maxITI, cuerewdelay, 1, postrewdelay);
    [eventlog_noinh] = simulateEvents(n_cues, 1, 2, 1, nan, ...
                meanITI, maxITI, cuerewdelay, 1, postrewdelay);
    % Inhibit during reward
    [optolog, first_inhib] = simulateInhibitionBlock(eventlog, 2, nan, ...
        -0.5, 1, inh_start, nan);
    
    % Calculate model values
    [DA,ANCCR,PRC, SRC, NC, Rs] = calculateANCCR(eventlog, T, ...
        alpha, k, samplingperiod, w,threshold,minimumrate,beta, ...
        alpha_r, maximumjitter, optolog, nan);

    % Calculate model values
    [DA_noinh,ANCCR_noinh,PRC_noinh, SRC_noinh, NC_noinh, Rs_noinh] = calculateANCCR(eventlog_noinh, T, ...
        alpha, k, samplingperiod, w,threshold,minimumrate,beta, ...
        alpha_r, maximumjitter, nan, nan);


    q_src = SRC.*Rs;
    q_src_noinh = SRC_noinh.*Rs_noinh;
    c1_idxs = find(eventlog(:,1) == 1);
    c1_idxs_noinh = find(eventlog_noinh(:,1) == 1);

    SRC_all_iters(:, :, :, iIter) = SRC(:, :, c1_idxs);
    Rs_all_iters(:, :, :, iIter) = Rs(:, :, c1_idxs);
    q_src_all_iters(:,:,iIter) = [squeeze(q_src(1, 2, c1_idxs)) + cost, zeros(length(c1_idxs), 1)]';

    SRC_all_iters_noinh(:, :, :, iIter) = SRC_noinh(:, :, c1_idxs_noinh);
    Rs_all_iters_noinh(:, :, :, iIter) = Rs_noinh(:, :, c1_idxs_noinh);
    q_src_all_iters_noinh(:,:,iIter) = [squeeze(q_src_noinh(1, 2, c1_idxs_noinh)) + cost, zeros(length(c1_idxs_noinh), 1)]';
end

p_action_all = nan(n_cues, nIter);
p_action_all_noinh = nan(n_cues, nIter);

close all;
for iIter = 1:nIter
    p_actions = squeeze(softmax(1, q_src_all_iters(:,:,iIter), 0.2));
    p_action_all(:, iIter) = p_actions;
    p_actions_noinh = squeeze(softmax(1, q_src_all_iters_noinh(:,:,iIter), 0.2));
    p_action_all_noinh(:, iIter) = p_actions_noinh;
    %plot(p_actions, 'Color', [0.6, 0.6, 0.6], 'LineWidth', 0.35)
    %plot(p_actions_noinh, 'Color', [0.6 0.6 1], 'LineWidth', 0.35)
end

%% Plot

p_action_all_mean = mean(p_action_all');
x = 1:numel(p_action_all_mean);
sem = std(p_action_all') / sqrt(nIter);
curve1 = p_action_all_mean + sem;
curve2 = p_action_all_mean - sem;
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];

p_action_all_mean_noinh = mean(p_action_all_noinh');
x_noinh = 1:numel(p_action_all_mean_noinh);
sem_noinh = std(p_action_all_noinh') / sqrt(nIter);
curve1_noinh = p_action_all_mean_noinh + sem_noinh;
curve2_noinh = p_action_all_mean_noinh - sem_noinh;
x2_noinh = [x_noinh, fliplr(x_noinh)];
inBetween_noinh = [curve1_noinh, fliplr(curve2_noinh)];

dir = 'D:\OneDrive - University of California, San Francisco\dopamine contingency\figures\figures for r3';

fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 4 3.5]);
hold on;
fill(x2, inBetween, [1, 0, 0],'FaceAlpha',0.3,'EdgeColor','none');
fill(x2_noinh, inBetween_noinh', [0 0 0],'FaceAlpha',0.3,'EdgeColor','none');
plot(x, p_action_all_mean, 'Color', [1, 0, 0], 'LineWidth', 0.5)
plot(x_noinh, p_action_all_mean_noinh, 'Color', [0 0 0], 'LineWidth', 0.5)
xline(inh_start,'k:', 'LineWidth', 0.5)
xlabel('Trial')
ylabel('p(lick|cue)')

set(gca,'Box','off','TickDir','out','FontSize',8,'LineWidth',0.35,...
    'YTick',0:0.5:1,'YLim',[0 1],'XTick',0:200:400,'XLim',[0 400])

cd(dir);
print(fHandle, '-depsc', '-painters', 'extinction.ai')
close all

% Define softmax
function [pi] = softmax(i, q, T)
 % Adding 1 approximates a zero energy null action, only do this for
 % extinction
    pi = exp(q(i,:)/T) ./ sum(exp(q/T));
end



