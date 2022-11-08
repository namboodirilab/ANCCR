clearvars; clc; close all;

rng(2)

directory = 'D:\OneDrive - UCSF\Huijeong\DA';

mouseList = {'HJ_FP_M2';'HJ_FP_M3';'HJ_FP_M4';'HJ_FP_F1';'HJ_FP_F2';'HJ_FP_M6';'HJ_FP_M7'};
nMouse = length(mouseList);
startday = [16,12,13,7,13,7,6]; % first session
learnedday = [22,17,25,12,17,12,14]; % first session w/ anticipatory lick
endday = [28,23,31,18,23,18,19]; % last session
ndays = [endday-startday]+1;
[cuersp,day] = deal(cell(nMouse,1));
for iM = 1:nMouse
    iM
    % pool sessions for each animal
    behfile = findfiles(mouseList{iM},[directory,'\',mouseList{iM},'\Pavlovian'],1,'Day','.mat');
    days = cellfun(@(y) str2double(y(4)),cellfun(@(x) strsplit(fileparts(x),{'Day','_'}),...
        behfile,'UniformOutput',false));
    indays = days>=startday(iM) & days<=endday(iM);
    behfile = behfile(indays);
    [days,sortidx] = sort(days(indays));
    behfile = behfile(sortidx);
    day{iM} = days;

    cuersp{iM} = nan(length(behfile),1);
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
        [~,~,CS,CStime,~,~,fxreward] = eventfrompavlovian(eventtime,eventindex,nosolenoidflag);
        CSrw = unique(CS(fxreward==1));

        % calculate CS+ response
        [timecs,cssignal] = alignsignal2event(T(:,1),dff',CStime(CS==CSrw),[-2000 3000],0);
        cuersp{iM}(i) = mean(aucsignal(cssignal,timecs,[0 1000]));
        cuersp{iM} = rmmissing(cuersp{iM});
    end
end

%% FigS10C
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 5 3.5]);
learnedidx = cellfun(@(x,y) find(x==y),day,num2cell(learnedday)','UniformOutput',false);
data = cell2mat(cellfun(@(x,y) x([1,2,[-2:2]+y,end-1,end])/mean(x(end)),...
    cuersp,learnedidx,'UniformOutput',false));
plot(1:size(data,2),data','Color',[0.6 0.6 0.6],'LineWidth',0.35)
hold on
errorbar(1:size(data,2),mean(data),std(data)/sqrt(nMouse),'k')
plot([0 size(data,2)+1],[1 1],'k:','LineWidth',0.35);
ylabel({'Normalized CS response';'(1 = response of last session)'})
xlabel('Session')
set(gca,'Box','off','TickDir','out','FontSize',8,'LineWidth',0.35,...
    'XTick',1:size(data,2),'XTickLabel',{'1','2','x-2','x-1','x','x+1','x+2','n-1','n'},...
    'YLim',[-0.3 1.7],'YTick',[0:0.5:1.5],'XTickLabelRotation',45)
%%
mean(mean(data(:,[5]),2))

