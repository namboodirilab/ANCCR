clearvars; clc; close all;

directory = 'D:\OneDrive - UCSF\Huijeong\DA';

mouseList = {'HJ_FP_M2';'HJ_FP_M3';'HJ_FP_M4';'HJ_FP_F1';'HJ_FP_F2';'HJ_FP_M6';'HJ_FP_M7';'HJ_FP_M8'};
nMouse = length(mouseList);

[rcorr_iri,out_bgd,out_nolick] = deal(NaN(15,nMouse));
rcorr_trial = NaN(nMouse,1);
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

    [auc,rewardnum] = deal(cell(nFile,1));
    n = 0;
    for iF = 1:nFile
        load(behfile{iF});
        
        try
            load([fileparts(behfile{iF}),'\Photometry.mat'],'dff','T');
        catch
            auc{iF} = deal(NaN(length(bgdrwinterval),1));
            continue
        end

        % exclude rewards w/ <3s IRI
        out_bgd(iF,iM) = mean(bgdrwinterval<3000 & [1:length(bgdrwtime)]'>1);
        out_nolick(iF,iM) = mean(isnan(firstlicktime));
        out = bgdrwinterval<3000 | isnan(firstlicktime);
        out(1) = false;  

        auc_signal = aucsignal2event(T(:,1),dff(:),firstlicktime,[-500 1000],out);
        auc_baseline = aucsignal2event(T(:,1),dff(:),firstlicktime,[-2000 -500],out);
        % normalized dopamine response by subtracting baseline response
        auctemp =  auc_signal-auc_baseline;

        auc{iF} = auctemp(~out);
        rewardnum{iF} = find(~out)+n; 

        % correlation b/w dopamine response and iri  within each session
        rcorr_iri(iF,iM) = corr(bgdrwinterval(~out),auc{iF},'rows','complete');
        n = n+length(firstlicktime);
    end
    % correlation b/w dopamine response and reward number across sessions
    rcorr_trial(iM) = corr(cell2mat(rewardnum),cell2mat(auc),'rows','complete');
end

%% correlation b/w dopamine response and reward number 
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 1.65 2.9]);
hold on;
bar(0,nanmean(rcorr_trial),1,'FaceColor',[0.8 0.8 0.8],'EdgeColor','k','LineWidth',0.35);
errorbar(0,nanmean(rcorr_trial),nanstd(rcorr_trial)/sqrt(nMouse),'k','Linewidth',0.5);
scatter(rand(nMouse,1)*0.8-0.4,rcorr_trial,2.5,'k','filled');
ylim([0 0.8]);
xlim([-1 1]);
set(gca,'Box','off','TickDir','out','FontSize',7,'LineWidth',0.35,'XTick',[],'YTick',0:0.4:0.8);
ylabel('r (DA vs. # of rewards)','FontSize',8)
[~,pval,~,stat] = ttest(rcorr_trial);

%% correlation b/w dopamine response and iri 
r = nanmean(rcorr_iri,1); % session-average of each animal
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 1.65 2.9]);
hold on;
bar(0,nanmean(r),1,'FaceColor',[0.8 0.8 0.8],'EdgeColor','k','LineWidth',0.35);
errorbar(0,nanmean(r),nanstd(r)/sqrt(nMouse),'k','Linewidth',0.5);
scatter(rand(nMouse,1)*0.8-0.4,r,2.5,'k','filled');
ylim([0 0.6]);
xlim([-1 1]);
set(gca,'Box','off','TickDir','out','FontSize',7,'LineWidth',0.35,'XTick',[],'YTick',-0.2:0.2:0.6);
ylabel('r (DA vs. IRI)','FontSize',8)
[~,pval,~,stat] = ttest(r);

%% correlation b/w dopamine response and iri across sessions (figS8K)
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 4.5 3.5]);
hold on;
plot(1:15,rcorr_iri,'Color',[0.6 0.6 0.6],'LineWidth',0.35);
errorbar(1:7,nanmean(rcorr_iri(1:7,:),2),nanstd(rcorr_iri(1:7,:),[],2)./...
    sqrt(sum(~isnan(rcorr_iri(1:7,:)),2)),'k','LineWidth',0.5,'CapSize',3);
plot([0 12],[0 0],'k:');
ylim([-0.6 0.6]);
xlim([0 12])
set(gca,'Box','off','TickDir','out','FontSize',8,'LineWidth',0.35,'XTick',[1 5 9],'YTick',-0.6:0.3:0.6);
xlabel('Session');
ylabel('r (DA vs. previous IRI)');