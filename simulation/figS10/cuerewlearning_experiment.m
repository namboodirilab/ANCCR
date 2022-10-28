clearvars; clc; close all;
rng(2);

%% parameter set up
% task parameters
meanIRI = 12;
numrewards = 1000;

meanITI = 30;
cuerewdelay = [3 9; 9 9; 9 9]; % cue-duration change, background reward, extinction
outcomedelay = [3 3; 3 3; 3 3];
bgdrwperiod = [NaN NaN; NaN 6; NaN NaN];
mincuebgddelay = [0 0; 0 6; 0 0];
rwprob = [1 1; 1 1; 1 0];
numcue = 2000;
postrewdelay = 3;

% anccr model parameters
samplingperiod = 0.2;   
alpha_anccr = 0.02;        
alpha_r = 0.2;
w = 0.5;               
k = 1;                 
minimumrate = 10^(-3);
maximumjitter = 0.1;
beta = [0,0,1,0];
threshold = 0.6;
Tratio = 1.2;

% rpe model parameters - csc/microstimulus
alpha_rpe = 0.05;
gamma = 0.95;
lambda = 0;
statesize = 0.2;

nIter = 100;
%%
% row: CSC w/ ITI states, ANCCR, CSC w/o ITI states
% column: initial learning, cue duration change, background reward,
% extinction
cuersp = nan(3,size(cuerewdelay,1),nIter,round(numcue*4/3)); 
for iIter = 1:nIter
    iIter
    for iC = 1:size(cuerewdelay,1)
        eventlog_rr = simulateBackgroundRewards(numrewards,....
            meanIRI,3,1,0);
        [eventlog_c1,IRI_c1] = simulateEvents(repmat(numcue,1,2), [1,2], [3,nan], ...
            [1,nan], [4,nan], meanITI, meanITI*3, cuerewdelay(iC,1), [rwprob(iC,1),0],...
            postrewdelay, 3, bgdrwperiod(iC,1), mincuebgddelay(iC,1), 1);
        [eventlog_c2,IRI_c2] = simulateEvents(repmat(numcue,1,2), [1,2], [3,nan], ...
            [1,nan], [4,nan], meanITI, meanITI*3, cuerewdelay(iC,2), [rwprob(iC,2),0],...
            postrewdelay, 3, bgdrwperiod(iC,2), mincuebgddelay(iC,2), 1);
               
        for imdl = 1:2
            switch imdl
                case 1
                    eventlog = joinEventlogs(eventlog_c1,eventlog_c2);
                    [DA,~,eventtimeline] = simulateCSC(eventlog,3,statesize,alpha_rpe,gamma,lambda,[]);
                    incue = eventtimeline(:,1)==1;
                case 2
                    eventlog = joinEventlogs(eventlog_rr,eventlog_c1,eventlog_c2);
                    T = [ones(numrewards,1)*meanIRI; ones(size(eventlog_c1,1),1)*IRI_c1;...
                        ones(size(eventlog_c2,1),1)*IRI_c2];
                    T(isinf(T)) = max(T(~isinf(T)));
                    
                    [DA,ANCCR,PRC,SRC,NC] = calculateANCCR(eventlog, T*Tratio, alpha_anccr, k,...
                        samplingperiod,w,threshold,minimumrate,beta,alpha_r,maximumjitter,nan,[4,3],0);
                    incue = eventlog(:,1)==1;
            end
            cuersp_temp = DA(incue);
            
            if iC==1
                cuersp(imdl,1,iIter,1:numcue) = cuersp_temp(1:numcue); % initial learning
            end
            cuersp(imdl,iC+1,iIter,:) = cuersp_temp(end-round(numcue*4/3)+1:end);
        end
        
        if iC==2
            % w/o ITI states
            eventlog = joinEventlogs(eventlog_c1,eventlog_c2);
            [DA,~,eventtimeline] = simulateCSC(eventlog,3,statesize,alpha_rpe,gamma,lambda,[],0,cuerewdelay(iC,:));
            incue = eventtimeline(:,1)==1;
            cuersp_temp = DA(incue);
            cuersp(3,iC+1,iIter,:) = cuersp_temp(end-round(numcue*4/3)+1:end);
        end
    end
end

%%
cd('D:\OneDrive - University of California, San Francisco\figures\manuscript\dopamine_contingency\revision\data');
save('cuerewlearning_experiment.mat','cuersp');

%%
load('cuerewlearning_experiment.mat');
shortiri = load('cuerewlearning_shortiri.mat');
close all
x = [1,2,4,3];
test = [3,4,7,6];
mdl = {'RPE';'ANCCR'};
titlelist = {'initial learning','cue duration change','background reward','extinction'};
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 15 7]);
for imdl = 1:2
    for iC = 1:4
        axes('Position',axpt(4,2,x(iC),imdl,axpt(10,10,2:10,2:9)));
        data = rmmissing(squeeze(cuersp(imdl,iC,:,:)),2);
        data = cumsum(data,2)./repmat(abs(sum(data,2)),1,size(data,2));
        data_d = cell2mat(cellfun(@(x) decimate(x,10),mat2cell(data,...
            ones(size(data,1),1),size(data,2)),'UniformOutput',false));
        plot(decimate(1:size(data,2),10)/size(data,2),data_d,'Color',[0.8 0.8 0.8],'LineWidth',0.35);
        hold on;
        if iC==2 & imdl==2
            data2 = cumsum(shortiri.cuersp,2)./repmat(abs(sum(shortiri.cuersp,2)),1,size(shortiri.cuersp,2));
            data2_d = cell2mat(cellfun(@(x) decimate(x,10),mat2cell(data2,...
                ones(size(data2,1),1),size(data2,2)),'UniformOutput',false));
            plot(decimate(1:size(data2,2),10)/size(data2,2),data2_d,'Color',[1 0.8 0.8],'LineWidth',0.35);
            plot(decimate(1:size(data2,2),10)/size(data2,2),mean(data2_d),'Color','r','LineWidth',1);
        end
        plot(decimate(1:size(data,2),10)/size(data,2),mean(data_d),'Color','k','LineWidth',1);
        plot([0 1],[0 1],'k--','LineWidth',0.35)
        if iC>1
           plot(repmat(1/4,1,2),[-1 1.2],'k:','LineWidth',0.35);
        end
        if imdl~= 1 | iC~=3
           ylim([-0.1 1.1]) 
        else
            ylim([-1 1.1]);
            data = rmmissing(squeeze(cuersp(3,iC,:,:)),2);
            data = cumsum(data,2)./repmat(abs(sum(data,2)),1,size(data,2));
            data_d = cell2mat(cellfun(@(x) decimate(x,10),mat2cell(data,...
                ones(size(data,1),1),size(data,2)),'UniformOutput',false));
            plot(decimate(1:size(data,2),10)/size(data,2),data_d,'Color',[1 0.8 0.8],'LineWidth',0.35);
            hold on;
            plot(decimate(1:size(data,2),10)/size(data,2),mean(data_d),'Color','r','LineWidth',1);
        end
        
         set(gca,'Box','off','TickDir','out','FontSize',8,'XTick',0:0.5:1,'YTick',0:0.5:1);
         if iC>1
             set(gca,'YTickLabel',[]);
         else
             ylabel({mdl{imdl};['(Model ',num2str(imdl),')'];'Normalized cumsum (predicted DA)'});
         end
         
         if imdl==1
             set(gca,'XTickLabel',[]);
             title({['Experiment ',num2str(x(iC)+1)];['Test ',num2str(test(iC))];...
                 titlelist{iC}},'FontSize',8);
         else
             xlabel('Normalized trial');
         end
    end
end
%%
dir = 'D:\OneDrive - University of California, San Francisco\figures\manuscript\dopamine_contingency\revision';
print(fHandle,'-depsc','-painters',[dir,'\cuerewlearning_experiment.ai']);
