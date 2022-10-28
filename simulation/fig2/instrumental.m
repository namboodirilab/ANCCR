% Show extinction of DA response (also behavior)
clearvars; clc; 
rng(2)

%% Initialization
% Task parameters
meanITI = [30, 30];
maxITI = meanITI*3;
intercuedelay = [1, 1];
cuerewdelay = [3, 3];
postrewdelay = [2, 2];
n_cues = 10000;

% Model parameters
samplingperiod = 0.2;
alpha = 0.02;
alpha_r = 0.85;
w = 0.5;
k = 1;
Tratio = 1.2;
T = meanITI*Tratio;
minimumrate = 10^(-3);
threshold = 0.6;
maximumjitter = 0.1;
beta = [0; 0; 1];
Rtrue = [0; 0; 1];
exact_mean_or_not = 0;

% Inhibition parameters
inhib_targets = 3;
inhib_mcts = 1;
cue_labels = [1, 2];
rew_labels = [3, 3];
rew_mags = [1, 1];
r_probs = [0.9, 0.1];
inhib_mag = -1;
trials_off = 500;
trials_on = 1; % Unused

nIter = 20;
all_iter_p_action_c1_n = nan(16, nIter);
all_iter_p_action_c1_n1 = nan(16, nIter);

%% Simulation
% Simulate all cue and reward times
for iIter = 1:nIter
    [eventlog] = simulateEvents(n_cues*r_probs, cue_labels, rew_labels, ...
        rew_mags, [nan nan], meanITI, maxITI, cuerewdelay, r_probs, ...
        postrewdelay);
    
    % Inhibit during reward
    [optolog] = simulateInhibitionPattern(eventlog, ...
        inhib_targets, inhib_mcts, inhib_mag, trials_off);
    
    % Calculate model values
    [DA,ANCCR,PRC,SRC,NC, Rs] = calculateANCCR(eventlog, T(1), ...
        alpha, k, samplingperiod, w, threshold, minimumrate, beta, ...
        alpha_r, maximumjitter, optolog, nan);
    
    inh_r_idxs = find(optolog(:,1) == 1);
    inh_c1_idxs = find(optolog(:,1) == 2);
    c1_idxs = find(eventlog(:,1) == 1);
    
    da_c1_n = DA(inh_c1_idxs);
    c1_n1_idx = nan(length(inh_c1_idxs), 1);
    inh_count = 1;
    
    % Find n+1 idxs
    for i = 1:length(c1_idxs)
        if inh_count < length(inh_c1_idxs) + 1
            if c1_idxs(i) == inh_c1_idxs(inh_count)
                if i+1 <= length(c1_idxs)
                    c1_n1_idx(inh_count) = c1_idxs(i+1);
                    inh_count = inh_count + 1;
                else
                    % If there is no c1 that follows last incident of inh
                    c1_n1_idx(inh_count) = [];
                    inh_count = inh_count + 1;
                end
                
            end
        end
    end
    
    da_c1_n1 = DA(c1_n1_idx);
    
    q_src = SRC.*Rs;
    temperature = 0.5;
    p_action_c1_n = softmax(1,3,q_src(1:2,:,:), temperature, inh_c1_idxs);
    p_action_c1_n1 = softmax(1,3,q_src(1:2,:,:), temperature, c1_n1_idx);
    
    p_action_c1_n = rmmissing(squeeze(p_action_c1_n));
    p_action_c1_n1 = rmmissing(squeeze(p_action_c1_n1)); 

    all_iter_p_action_c1_n(:, iIter) = p_action_c1_n;
    all_iter_p_action_c1_n1(:, iIter) = p_action_c1_n1;
end
%% Plot
close all;
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 2.5 3.5]);
axes('Position',axpt(5,5,2:5,1:4)) 
hold on;
data = [mean(all_iter_p_action_c1_n,1)',mean(all_iter_p_action_c1_n1,1)'];
plot([1,2],data,'Color',[0.6 0.6 0.6], 'LineWidth',0.35)
errorbar([1,2], mean(data),std(data)/sqrt(nIter), 'Color', [0 0 0], "LineWidth",0.5)
set(gca,'Box','off','TickDir','out','FontSize',8,'LineWidth',0.35,...
    'XTick',1:2,'XTickLabel',{'n', 'n+1'},...
    'YTick',0:0.1:1,'YLim',[0.5 1],'XLim',[0.5 2.5])
ylabel('P(action 1)')
xlabel('Trial');

dir = 'D:\OneDrive - University of California, San Francisco\dopamine contingency\figures\figures for r3';
cd(dir);
print(fHandle,'-depsc','-painters','instrumental.ai')
% save('instrumental.mat', 'all_iter_p_action_c1_n', 'all_iter_p_action_c1_n1')

all_iter_p_action_c1_n_flat = mean(all_iter_p_action_c1_n, 1);
all_iter_p_action_c1_n1_flat = mean(all_iter_p_action_c1_n1, 1);
[significant_binary, p_value, ~, stats] = ttest(all_iter_p_action_c1_n_flat, all_iter_p_action_c1_n1_flat);

% Define softmax
function [pi] = softmax(i, j, q, T, idxs)
    pi = exp(q(i,j,idxs)/T) ./ sum(exp(q(:,j,idxs)/T));
end