clearvars; clc; close all;

rng(2)

directory = 'D:\OneDrive - UCSF\Huijeong\DA';

mouseList = {'HJ_FP_M2';'HJ_FP_M3';'HJ_FP_M4';'HJ_FP_F1';'HJ_FP_F2';'HJ_FP_M6';'HJ_FP_M7'};
nMouse = length(mouseList);

startday = [16,12,13,7,13,7,6]; % first session
endday = [28,23,31,18,23,18,19]; % last session
ndays = [endday-startday]+1;

psth = cell(nMouse,1);
[early_auc,late_auc] = deal(nan(nMouse,50*max(ndays))); % nmouse, ntrials (50*nsessions)
rw_auc = nan(nMouse,50);

for iM = 1:nMouse
    iM
    % pool sessions for each animal
    behfile = findfiles(mouseList{iM},[directory,'\',mouseList{iM},'\Pavlovian'],1,'Day','.mat');
    days = cellfun(@(y) str2double(y(4)),cellfun(@(x) strsplit(fileparts(x),{'Day','_'}),...
        behfile,'UniformOutput',false));
    indays = days>=startday(iM) & days<=endday(iM);
    behfile = behfile(indays);
    [~,sortidx] = sort(days(indays));
    behfile = behfile(sortidx);
    
    for i = 1:length(behfile)
        load(behfile{i});
        
        try
            load([fileparts(behfile{i}),'\Photometry.mat']);
        catch
            continue;
        end

        % load task related information
        eventtime = eventlog(:,2);
        eventindex = eventlog(:,1);
        nosolenoidflag = eventlog(:,3);
        [~,~,CS,CStime,~,~,fxreward,~,firstfxlicktime] =...
            eventfrompavlovian(eventtime,eventindex,nosolenoidflag);
        CSrw = unique(CS(fxreward==1));   

        % align dff to CS+ onset
        [timecs,cssignal] = alignsignal2event(T(:,1),dff',CStime(CS==CSrw),[-2000 3000],0);
        psth{iM} = [psth{iM};cssignal];
        
        if i==1
            % in first session, calculate reward response
            [timerw,rwsignal] = alignsignal2event(T(:,1),dff',firstfxlicktime(CS==CSrw),[-2000 3000],0);
            rw_auc(iM,1:sum(CS==CSrw)) = aucsignal(rwsignal,timerw,[0 1000]);
        end
    end
    % calculate AUC during early (0-1s) and late (2-3s) from cue onset
    early_auc(iM,1:size(psth{iM},1)) = aucsignal(psth{iM},timecs,[0 1000]); 
    late_auc(iM,1:size(psth{iM},1)) = aucsignal(psth{iM},timecs,[2000 3000]); 
end

%% Fig6B left
close all
dir = 'D:\OneDrive - University of California, San Francisco\figures\manuscript\dopamine_contingency\revision\fig6_new';
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 4 3]);
axes('Position',axpt(5,5,2:5,1:4));
plot(1:400,early_auc(:,1:400)./repmat(nanmean(early_auc(:,[-49:0]+400),2),1,400),'Color',[0.6 0.6 0.6],'LineWidth',0.35);
hold on;
plot(1:400,late_auc(:,1:400)./repmat(nanmean(early_auc(:,[-49:0]+400),2),1,400),'Color',[1 0.6 0.6],'LineWidth',0.35);
plot(1:400,mean(early_auc(:,1:400)./repmat(nanmean(early_auc(:,[-49:0]+400),2),1,400)),'Color','k','LineWidth',1);
plot(1:400,mean(late_auc(:,1:400)./repmat(nanmean(early_auc(:,[-49:0]+400),2),1,400)),'Color','r','LineWidth',1);
xlabel('Trial');
ylabel('Norm. DA response');
plot([10 30],[2 2],'Color','k','LineWidth',1);
plot([10 30],[1.8 1.8],'Color','r','LineWidth',1);
text(40,2,'Early (0-1 s)');
text(40,1.8,'Late (2-3 s)');
set(gca,'FontSize',8,'Box','off','TickDir','out','LineWidth',0.35,'XTick',0:100:400,...
    'YTick',-1:2,'YLim',[-1 2.3],'XTickLabelRotation',45)

%% Fig6B right
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 2.25 3]);
axes('Position',axpt(5,10,4:5,1:9))
hold on;
delta_auc = mean(late_auc(:,1:100)-early_auc(:,1:100),2)./nanmean(early_auc(:,[-49:0]+400),2);
bar(0,mean(delta_auc),1,'FaceColor',[0.6 0.6 0.6]);
errorbar(0,mean(delta_auc),std(delta_auc)/sqrt(nMouse),'k');
scatter(rand(nMouse,1)*0.8-0.4,delta_auc,2,'k','Filled');
ylabel({'\DeltaNorm. DA response';'(late-early)'});
set(gca,'FontSize',8,'XTickLabel',[],'XLim',[-1 1],'Box','off','TickDir','out',...
    'LineWidth',0.35,'YTick',-0.2:0.1:0.2,'YLim',[-0.2 0.2])
[~,p,~,stat] = ttest(delta_auc)
