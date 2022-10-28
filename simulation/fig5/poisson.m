clearvars; clc; close all

rng(2)

meanITI = [30 33];
maxITI = meanITI*3;
cuerewdelay = [9 3];
outcomedelay = 3;
cueduration = 0.25;
numcue = [2000, 1000];

samplingperiod = 0.2;
alpha = 0.02;
alpha_r = 0.2;
w = 0.5;
k = 1;
threshold = 0.6;
beta = [0,0,1]';
maximumjitter = 0.1;
minimumrate = 10^(-3);
Tratio = 1.2;

% rpe model parameters - csc/microstimulus
alpha_rpe = 0.05;
gamma = 0.95;
lambda = [0, 0.95];
statesize = 0.2;
sigma = 0.08;
nmicrostimulus = 20;
d = 0.99;
nIter = 100;


%% Simulate experiment
[cueave,rwave] = deal(nan(nIter,2,3));
for iiter = 1:nIter
    iiter
    
    % pre-conditioning w/ trial structure
    eventlog_cond = simulateEvents(repmat(numcue(1),1,2),[1,2],3,...
        [1,1],[nan,nan],repmat(meanITI(1),1,2),repmat(maxITI(1),1,2),...
        repmat(cuerewdelay(1),1,2),[1,0],repmat(outcomedelay,1,2));
    
    % trial-less conditioning
    eventlog = simulateEventsTrialLess(numcue(2), 1, 3, 1, nan,...
        meanITI(2), maxITI(2), cueduration, cuerewdelay(2), 1);
    eventlog(:,2) = eventlog(:,2)+eventlog_cond(end,2);
    eventlog = [eventlog_cond; eventlog];
    
    T = meanITI(2)*Tratio;

    for imdl = 1:3
        switch imdl
            case 1
                [DA,~,eventtimeline] = simulateCSC(eventlog,3,statesize,alpha_rpe,gamma,lambda(1),[]);
            case 2
                [DA,~,eventtimeline] = simulatemicrostimulus(eventlog,3,nmicrostimulus,statesize,...
                    alpha_rpe,gamma,lambda(2),sigma,d,[]);
            case 3
                DA = calculateANCCR(eventlog, T, alpha, k,...
                    samplingperiod,w,threshold,minimumrate,beta,alpha_r,maximumjitter,nan,nan);
        end
        
        if imdl<3
            incue = find(eventtimeline(:,1)==1,numcue(2),'last');
            cuetimes = incue*statesize;
            cuersp = DA(incue);
            rwrsp = DA(incue+round(cuerewdelay(2)/statesize));
        else
            incue = find(eventlog(:,1)==1,numcue(2),'last');
            cuetimes = eventlog(incue,2);
            cuersp = DA(incue);
            rwrsp = DA(find(eventlog(:,1)==3,numcue(2),'last'));
        end
        
        % finding intermediate cues and previous cues of them
        % if delay b/w pair is shorter than 0.5, excluded it from analysis
        intermediateidx = [false;diff(cuetimes)<cuerewdelay(2) & diff(cuetimes)>=0.5];
        trialidx = [[intermediateidx(2:end);0],intermediateidx]; 

        % if the delays from previous cue and to next cue are both shorter
        % than cuerewdelay, excluded all three cues from analysis
        trialidx(find(sum(trialidx,2)==2)-1,:) = 0;
        trialidx(find(sum(trialidx,2)==2)+1,:) = 0;
        trialidx(sum(trialidx,2)==2,:) = 0;
        
        for i = 1:2
            cueave(iiter,i,imdl) = mean(cuersp(trialidx(:,i)==1));
            rwave(iiter,i,imdl) = mean(rwrsp(trialidx(:,i)==1));
        end        
    end
end

%%
cd('D:\OneDrive - University of California, San Francisco\figures\manuscript\dopamine_contingency\revision\data');
save('poisson.mat','cueave','rwave');

%%
x = [0.5 2 4.5];
typelist = {'cue';'reward'};
model = {'RPE (CSC)','RPE (MS)','ANCCR'};
for i = 1:2
    if i==1
        ratio = squeeze(cueave(:,2,:)./cueave(:,1,:));
    else
        ratio = squeeze(rwave(:,2,:)./rwave(:,1,:));
    end
    fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 3.7 4.5]);
    hold on;
    
    for ii = 1:3
    bar(x(ii),nanmean(ratio(:,ii)),1,'FaceColor',[0.8 0.8 0.8],'EdgeColor','k','LineWidth',0.35);
    [~,p,~,stat] = ttest(ratio(:,ii),i-1)
    end
    errorbar([0.5 2 4.5],nanmean(ratio),nanstd(ratio)/sqrt(nIter),...
        'k','Linewidth',0.5,'CapSize',3,'LineStyle','none');
  
    set(gca,'XTick',[0.5 2 4.5],'XTickLabel',model,'XTickLabelRotation',45,...
        'Box','off','TickDir','out','FontSize',8,'LineWidth',0.35);
    ylabel({'Ratio of predicted';[typelist{i},' response'];'(intermediate / previous)'});
    if i==1
        ylim([-1 1]);
        set(gca,'YTick',-1:1);
    else
        plot([-0.5 5.5],[1 1],'k:','LineWidth',0.35);
        ylim([0 2.5]);
        set(gca,'YTick',0:2);
    end
    xlim([-0.5 5.5]);
        cd('D:\OneDrive - University of California, San Francisco\figures\manuscript\dopamine_contingency\fig4');
        print(fHandle,'-depsc',['model_',typelist{i},'.ai']);
end