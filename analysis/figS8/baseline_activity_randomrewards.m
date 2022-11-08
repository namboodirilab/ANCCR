clearvars; clc; close all;

directory = 'D:\OneDrive - UCSF\Huijeong\DA';

mouseList = {'HJ_FP_M2';'HJ_FP_M3';'HJ_FP_M4';'HJ_FP_F1';'HJ_FP_F2';'HJ_FP_M6';'HJ_FP_M7';'HJ_FP_M8'};
nMouse = length(mouseList);

% event indices
index_lick = 5;
index_bgdrw = 7;
index_sessionend = 0;

rcorr = nan(nMouse,2);
for iM = 1:nMouse
    iM
    % pool behavior files
    behfile = findfiles(mouseList{iM},[directory,'\',mouseList{iM},'\Randomrewards'],1,'Day','.mat');
    days = cellfun(@(y) str2double(y(4)),cellfun(@(x) strsplit(fileparts(x),{'Day','_'}),...
        behfile,'UniformOutput',false));
    [~,sortidx] = sort(days);
    behfile = behfile(sortidx);
    index_30s = cellfun(@(x) contains(fileparts(x),'30s'),behfile);
    behfile(find(index_30s,1,'first'):end) = []; % use 12s IRI sessions
    nFile = length(behfile);

    [aucrw,baselinepk] = deal(NaN(nFile,1));
    for iF = 1:nFile
        load(behfile{iF});
        try
            load([fileparts(behfile{iF}),'\Photometry.mat'],'dff','T');
        catch
            continue
        end

        eventtime = eventlog(:,2);
        eventindex = eventlog(:,1); 

        % calculate event time stamps
        licktime = eventtime(eventindex==index_lick);
        bgdrwtime = eventtime(eventindex==index_bgdrw);
        bgdrwinterval = [NaN;diff(bgdrwtime)];
        sessionendtime = eventtime(eventindex==index_sessionend);
        firstlicktime = firsttimeafterevent(licktime,bgdrwtime,[bgdrwtime(2:end);sessionendtime]);

        % set threshold in the first session of each animal
        if iF==1
            pks = findpeaks(dff*100,'MinPeakDistance',50);
            sorted = sort(pks,'descend');
            threshold = nanmean(sorted(1:length(sorted)*0.05))*0.2;
        end
        
        % calculate reward response
        shortiri = bgdrwinterval<3000; % exclude rewards w/ <3s iri
        shortiri(1) = false;
        % pre-reward baseline activity
        aucbase = aucsignal2event(T(:,1),dff,firstlicktime(~shortiri),[-2000 -500]); 
        % baseline-corrected reward response
        aucrw(iF) =  nanmean(aucsignal2event(T(:,1),dff,firstlicktime(~shortiri),[-500 1000])-aucbase); 
        

        % find baseline peaks
        [pks,locs] = findpeaks(dff*100,'MinPeakHeight',threshold,...
            'MinPeakDistance',50,'MinPeakProminence',threshold*0.5);
        pktime = T(locs,1);
       
        % exclude peaks during reward consumption ([0 2] s from firstlicktime) 
        infirstlick = sum(cell2mat(cellfun(@(x) pktime>=x & pktime<x+2000,...
            num2cell(firstlicktime),'UniformOutput',false)'),2);
        
        % exclude peaks w/ <0.5s latency from any lick
        latencyfromlastlick = cellfun(@(x) x-licktime(find(x-licktime>0,1,'last')),num2cell(pktime),'UniformOutput',false);
        latencytonextlick =  cellfun(@(x) licktime(find(x-licktime<0,1,'first'))-x,num2cell(pktime),'UniformOutput',false);
        latencyfromlastlick(cellfun(@isempty,latencyfromlastlick)) = {NaN};
        latencytonextlick(cellfun(@isempty,latencytonextlick)) = {NaN};
        aroundlick = cell2mat(latencyfromlastlick)<500 | cell2mat(latencytonextlick)<500;
        
        baselinepk(iF) = mean(pks(~infirstlick & ~aroundlick));
    end
    rcorr(iM,1) = corr([1:nFile]',aucrw); % correlation b/w session number and reward response
    rcorr(iM,2) = corr([1:nFile]',baselinepk); % correlation b/w session number and basline peak size
end

%% correlation b/w DA baseline activity vs. session number (FigS8G)
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 2.2 3.5]);
hold on;
bar(0,mean(rcorr(:,2)),1,'FaceColor',[0.8 0.8 0.8],'EdgeColor','k','LineWidth',0.35);
errorbar(0,mean(rcorr(:,2)),std(rcorr(:,2))/sqrt(nMouse),'k','Linewidth',0.5);
scatter(rand(nMouse,1)*0.8-0.4,rcorr(:,2),2.5,'k','filled');
ylabel({'r (DA baseline activity'; 'vs. session)'});
set(gca,'Box','off','TickDir','out','FontSize',8,'LineWidth',0.35,...
    'XTick',[],'YTick',[-1 0 1],'XLim',[-1 1]);

%% correlation b/w DA baseline activity vs. reward response (FigS8H)
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
