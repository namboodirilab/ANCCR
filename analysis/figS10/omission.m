clearvars; clc; 

rng(7)

directory = 'D:\OneDrive - UCSF\Huijeong\DA';

mouseList = {'HJ_FP_M2';'HJ_FP_M3';'HJ_FP_M4';'HJ_FP_F1';'HJ_FP_F2';'HJ_FP_M6';'HJ_FP_M7'};
nMouse = length(mouseList);

% the first day of each condition; cue duration change and extinction
startday = [29,24,32,19,24,19,20; 70,61,66,42,46,47,43];

aucdata = cell(nMouse,2);
% analysis window expecting to see omission response: 3-4 s from cue onset
% in cue duration change, 9-10 s from cue onset in extinction 
win = [3000 4000; 9000 10000];
for iM = 1:nMouse
    iM
    behfile = findfiles(mouseList{iM},[directory,'\',mouseList{iM},'\Pavlovian'],1,'Day','.mat');
    days = cellfun(@(y) str2double(y(4)),cellfun(@(x) strsplit(fileparts(x),{'Day','_'}),...
        behfile,'UniformOutput',false));
    for i = 1:2
        in = days==startday(i,iM);
        load(behfile{in})

        % load task related information
        eventtime = eventlog(:,2);
        eventindex = eventlog(:,1);
        nosolenoidflag = eventlog(:,3);
        [~,~,CS,CStime,~,~,fxreward] = eventfrompavlovian(eventtime,eventindex,nosolenoidflag);
        if i==1
            CSrw = unique(CS(fxreward==1));
        end

        load([fileparts(behfile{in}),'\Photometry.mat']);
        
        % calculate AUC during analysis window
        [timecs,cssignal] = alignsignal2event(T(:,1),dff',CStime(CS==CSrw),[-2000 12000],10);
        aucbase = aucsignal(cssignal,timecs,[-1000 0]);
        auccs = aucsignal(cssignal,timecs,win(i,:))/diff(win(i,:));

        % normalize response by pre-cue baseline response
        aucdata{iM,i} = auccs-aucbase;
    end
end

%% FigS10B
aucdata_cumsum = cellfun(@cumsum,aucdata,'UniformOutput',false);
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 4 3.5]);
hold on;
for iM = 1:nMouse
    plot(aucdata_cumsum{iM,1},'k')
    hold on;
    plot(aucdata_cumsum{iM,2},'r');
end
[~,p,~,stat] = ttest(cellfun(@(x) x(end),aucdata_cumsum(:,1)),cellfun(@(x) x(end),aucdata_cumsum(:,2)));
xlabel('Trial');
ylabel('Cumsum (DA response)');
set(gca,'Box','off','TickDir','out','FontSize',8,'XTick',[0 25 50],'YTick',-200:100:200);
xlim([0 51])
ylim([-200 200]);

