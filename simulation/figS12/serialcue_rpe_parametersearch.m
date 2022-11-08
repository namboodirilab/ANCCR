% simulate sequential learning task with different lambda

clearvars; close all; clc;
rng(2)

% task parameter
numcue = 800;
cuerewdelay = 3;
cuecuedelay = 1.5;
consumdelay = 3;
meanITI = 60;

inhibitionwrtcue = [cuecuedelay, cuerewdelay];

% model parameter
alpha = 0.05; % learning rate
gamma = 0.95; % discount factor
lambda = [0,0.5,0.95];
numstimulus = 20; % number of microstimulus
d = 0.99; % memory decay constant
statesize = 0.2;
sigma = 0.08; % width of basis function

nExp = size(meanITI,1); % number of experiments

%%

nIter = 20;
darsp = cell(2,length(lambda),3);
valuemap = cell(2,length(lambda));
for iIter = 1:nIter
    % generate eventlog
    eventlog = simulateEventChain(numcue, 2, 4, meanITI, ...
        meanITI*3, cuecuedelay, cuerewdelay, 1, consumdelay);
    
    for irpe = 1:2 % 1: csc, 2: microstimulus
        for iparam = 1:length(lambda) % different lambda
            printf('%d iteration: %d lambad, %d rpe model\n',iIter,iparam,irpe);
            if irpe==1
                [rpetimeline,valuetimeline,eventtimeline,inhibitiontimeline] =...
                    simulateCSC(eventlog,3,statesize,alpha,gamma,lambda(iparam),[]);
            else
                [rpetimeline,valuetimeline,eventtimeline,inhibitiontimeline] =...
                    simulatemicrostimulus(eventlog,3,numstimulus,...
                    statesize,alpha,gamma,lambda(iparam),sigma,d,[]);
            end
            
            % pool DA response
            for ie = 1:3 % 1:cs1, 2:cs2, 3:reward
                inevent = eventtimeline(:,1)==ie;
                darsp{irpe,iparam,ie}(:,iIter) = rpetimeline(inevent);
            end
        end
    end
end
%%

%% FigS12A
close all
modelList = {'CSC';'Microstimulus'};
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 12 6]);
for irpe = 1:2
    for iparam = 1:length(lambda)
        axes('Position',axpt(length(lambda),2,iparam,irpe))
        hold on;
        plot(darsp{irpe,iparam,1},'Color',[0.6 0.6 0.6],'LineWidth',0.35)
        plot(darsp{irpe,iparam,2},'Color',[0.6 0.6 1],'LineWidth',0.35)

        plot(mean(darsp{irpe,iparam,1},2),'k','LineWidth',1)
        plot(mean(darsp{irpe,iparam,2},2),'b','LineWidth',1)

        ylim([-0.15 0.4])
        set(gca,'Box','off','TickDir','out','FontSize',8,'LineWidth',0.35,...
            'XTick',0:200:800,'XTickLabelRotation',45,'YTick',-0.2:0.2:0.4);
        if irpe==1
            title(['Lambda = ',num2str(lambda(iparam))])
            set(gca,'XTickLabel',[]);
        else
            xlabel('Trial')
        end
        if iparam==1
            ylabel({modelList{irpe};'Predicted DA response'})
        else
            set(gca,'YTickLabel',[])
        end
        
    end
end

%% FigS12B
close all
aveda = cellfun(@(x) mean(x(end-49:end,:)),darsp(:,:,:),'UniformOutput',false);
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 4 7]);
for irpe = 1:2
    axes('Position',axpt(1,2,1,irpe,axpt(5,5,3:5,1:4)))
    hold on;
    for iparam = 1:length(lambda)
        scatter(rand(nIter,1)+(iparam-1)*1.5,aveda{irpe,iparam,1},2,[0.6 0.6 0.6],'filled')
        scatter(rand(nIter,1)+(iparam-1)*1.5,aveda{irpe,iparam,2},2,[0.6 0.6 1],'filled')
        errorbar((iparam-1)*1.5+0.5,mean(aveda{irpe,iparam,1}),std(aveda{irpe,iparam,1})/sqrt(nIter),'Color','k','LineWidth',0.5)
        errorbar((iparam-1)*1.5+0.5,mean(aveda{irpe,iparam,2}),std(aveda{irpe,iparam,2})/sqrt(nIter),'Color','b','LineWidth',0.5)
    end
    plot([-1 5],[0 0],'k:')
    ylabel({modelList{irpe};'Predicted DA response';'during last 50 trials'})
    xlim([-1 5])
    ylim([-0.05 0.4])
    set(gca,'Box','off','TickDir','out','FontSize',8,'LineWidth',0.35,'XTick',0.5:1.5:3.5)
    if irpe==1
        set(gca,'XTickLabel',[]);
    else
        set(gca,'XTickLabel',lambda)
        xlabel('Lambda')
    end
end
print(fHandle,'-depsc','-painters','D:\OneDrive - UCSF\figures\manuscript\dopamine_contingency\revision\figS12_new\rpe_simulation_lambdas_last50.ai')

%%
