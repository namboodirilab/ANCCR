clearvars; clc; close all;

directory = 'D:\OneDrive - UCSF\Huijeong\DA';

mouseList = {'HJ_FP_M2';'HJ_FP_M3';'HJ_FP_M4';'HJ_FP_F1';'HJ_FP_F2';'HJ_FP_M6';'HJ_FP_M7';'HJ_FP_M8'};
nMouse = length(mouseList);

rcorr = nan(nMouse,2);
for iM = 1:nMouse
    iM
    behfile = findfiles('Events_randomrewards.mat',[directory,'\',mouseList{iM},'\Randomrewards'],1,'Day');
    days = cellfun(@(y) str2double(y(4)),cellfun(@(x) strsplit(fileparts(x),{'Day','_'}),...
        behfile,'UniformOutput',false));
    [~,sortidx] = sort(days);
    behfile = behfile(sortidx);
    index_30s = cellfun(@(x) contains(fileparts(x),'30s'),behfile);
    behfile(find(index_30s,1,'first'):end) = []; % use only 12s IRI sessions
    nFile = length(behfile);

    [aucrw,baselinepk] = deal(NaN(nFile,1));
    for iF = 1:nFile
        load(behfile{iF});
        try
            load([fileparts(behfile{iF}),'\Photometry.mat'],'dff','T');
        catch
            continue
        end
        shortiri = bgdrwinterval<3000;
        shortiri(1) = false;

        % set threshold in the first session of each animal
        if iF==1
            pks = findpeaks(dff*100,'MinPeakDistance',50);
            sorted = sort(pks,'descend');
            threshold = nanmean(sorted(1:length(sorted)*0.05))*0.2;
        end
        
        % reward response
        aucbase = aucsignal2event(T(:,1),dff,firstlicktime(~shortiri),[-2000 -500]);
        aucrw(iF) =  nanmean(aucsignal2event(T(:,1),dff,firstlicktime(~shortiri),[-500 1000])-aucbase);
        

        % baseline peaks
        [pks,locs] = findpeaks(dff*100,'MinPeakHeight',threshold,...
            'MinPeakDistance',50,'MinPeakProminence',threshold*0.5);
        pktime = T(locs,1);
       
        infirstlick = sum(cell2mat(cellfun(@(x) pktime>=x & pktime<x+2000,...
            num2cell(firstlicktime),'UniformOutput',false)'),2);
        
        latencyfromlastlick = cellfun(@(x) x-licktime(find(x-licktime>0,1,'last')),num2cell(pktime),'UniformOutput',false);
        latencytonextlick =  cellfun(@(x) licktime(find(x-licktime<0,1,'first'))-x,num2cell(pktime),'UniformOutput',false);
        latencyfromlastlick(cellfun(@isempty,latencyfromlastlick)) = {NaN};
        latencytonextlick(cellfun(@isempty,latencytonextlick)) = {NaN};
        aroundlick = cell2mat(latencyfromlastlick)<500 | cell2mat(latencytonextlick)<500;
        
        baselinepk(iF) = mean(pks(~infirstlick & ~aroundlick));
    end
    rcorr(iM,1) = corr([1:nFile]',aucrw);
    rcorr(iM,2) = corr([1:nFile]',baselinepk);
end

%% correlation b/w DA baseline activity vs. session number
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 2.2 3.5]);
hold on;
bar(0,mean(rcorr(:,2)),1,'FaceColor',[0.8 0.8 0.8],'EdgeColor','k','LineWidth',0.35);
errorbar(0,mean(rcorr(:,2)),std(rcorr(:,2))/sqrt(nMouse),'k','Linewidth',0.5);
scatter(rand(nMouse,1)*0.8-0.4,rcorr(:,2),2.5,'k','filled');
ylabel({'r (DA baseline activity'; 'vs. session)'});
set(gca,'Box','off','TickDir','out','FontSize',8,'LineWidth',0.35,...
    'XTick',[],'YTick',[-1 0 1],'XLim',[-1 1]);

%% correlation b/w DA baseline activity vs. reward response
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 4 3.8]);
hold on;
scatter(rcorr(:,2),rcorr(:,1),3.5,'k','filled');
beta = glmfit(rcorr(:,2),rcorr(:,1));
plot([-1 1],[-1 1]*beta(2)+beta(1),'k');
plot([-1 1 NaN 0 0],[0 0 NaN -1 1],'k:');
[rr,pp] = corr(rcorr(:,2),rcorr(:,1));
text(-0.9,-0.8,['r = ',num2str(round(rr*100)/100)],'FontSize',7);
text(-0.9,-0.9,['p = ',num2str(round(pp*1000)/1000)],'FontSize',7);
xlabel({'r (DA baseline activity'; 'vs. session)'});
ylabel({'r (DA reward response'; 'vs. session)'});
set(gca,'Box','off','TickDir','out','FontSize',8,'LineWidth',0.35,...
    'XTick',[-1 0 1],'YTick',[-1 0 1]);
