clearvars; clc; close all;

directory = 'D:\OneDrive - UCSF\Huijeong\DA';

mouseList = {'HJ_FP_M2';'HJ_FP_M3';'HJ_FP_M4';'HJ_FP_F1';'HJ_FP_F2';'HJ_FP_M6';'HJ_FP_M7';'HJ_FP_M8'};
nMouse = length(mouseList);

lickrate = nan(nMouse,15,3);
rcorr = nan(nMouse,1);
[avepsth_nonconsum,avepsth_consum] = deal(nan(nMouse,15,1085));
for iM = 1:nMouse

    behfile = findfiles('Events_randomrewards.mat',[directory,'\',mouseList{iM},'\Randomrewards'],1,'Day');
    days = cellfun(@(y) str2double(y(4)),cellfun(@(x) strsplit(fileparts(x),{'Day','_'}),...
        behfile,'UniformOutput',false));
    [days,sortidx] = sort(days);
    behfile = behfile(sortidx);
    index_30s = cellfun(@(x) contains(fileparts(x),'30s'),behfile);
    behfile(find(index_30s,1,'first'):end) = []; % use only 12s IRI sessions
    nFile = length(behfile);
    
    [aucrw_total,lickrw_total] = deal([]);
    for iF = 1:nFile
        load(behfile{iF});
        try
            load([fileparts(behfile{iF}),'\Photometry.mat'],'dff','T');
        catch
            aucrw_total = [aucrw_total;NaN(length(bgdrwtime),1)];
            continue
        end

        %% averaged lick rates of session
        lickrate(iM,iF,1) = length(lickconsumidx)/(sessionendtime/1000); % total lick
        lickrate(iM,iF,2) = sum(lickconsumidx)/...
            (sum(consumboutlength(~isnan(consumboutlength)))/1000); % consummatory lick
        lickrate(iM,iF,3) = sum(~lickconsumidx)/...
            ((sessionendtime-sum(consumboutlength(~isnan(consumboutlength))))/1000); %non-consummatory lick

        %% reward-by-reward lick rate and response
        % exclued rewards if interval from previous reward is less than 3 s
        % or lick bout was continued from previous reward
        shortiri = bgdrwinterval<3000;
        shortiri(1) = false;
        separatebout = cellfun(@(x) sum(ismember([1,2],lickboutidx(lickrwidx==x)))==2,...
            num2cell(1:length(bgdrwtime)))';

        aucbase = aucsignal2event(T(:,1),dff,firstlicktime,[-2000 -500]);
        aucrw = aucsignal2event(T(:,1),dff,firstlicktime,[-500 1000])-aucbase;
        aucrw_total = [aucrw_total;aucrw(~shortiri & separatebout)];

        lickrw = 1000*cellfun(@(x,y) sum(licktime>=x & licktime<y & lickconsumidx),...
            num2cell(bgdrwtime),num2cell([bgdrwtime(2:end);sessionendtime]))./consumboutlength;
        lickrw_total = [lickrw_total; lickrw(~shortiri & separatebout)];

        %% averaged psth of non-consummatory vs. consummatory lick
        nonconsum_lick = licktime(lickboutidx==1 & ~lickconsumidx);
        [time,psth_nonconsum] =...
            alignsignal2event(T(:,1),dff(:)*1000,nonconsum_lick,[-3 6]*10^3,0);
        [~,psth_consum] =...
            alignsignal2event(T(:,1),dff(:)*1000,firstlicktime,[-3 6]*10^3,0);

        psth_baseline = nanmean(psth_nonconsum(:,time>=-1500 & time<-500),2);
        avepsth_nonconsum(iM,iF,:) = (nanmean(psth_nonconsum)-nanmean(psth_baseline))/nanstd(psth_baseline);
        
        psth_baseline = nanmean(psth_consum(:,time>=-1500 & time<-500),2);
        avepsth_consum(iM,iF,:) = (nanmean(psth_consum)-nanmean(psth_baseline))/nanstd(psth_baseline);
    end
    out = isinf(lickrw_total);
    rcorr(iM) = corr(lickrw_total(~out),aucrw_total(~out),'rows','complete');
end

%% lick rate acorss sessions
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 7 4]);
titleList = {'total';'consumption';'non-consumption'};
for i = 1:3
h(i) = subplot(1,3,i); hold on;
plot(1:15,squeeze(lickrate(:,:,i)),'Color',[0.6 0.6 0.6]);
errorbar(1:7,nanmean(lickrate(:,1:7,i),1),nanstd(lickrate(:,1:7,i),[],1)./...
    sqrt(sum(~isnan(lickrate(:,1:7,i)),1)),'k','Capsize',3);
title(h(i),titleList{i});
end
set(h,'YLim',[0 8],'XLim',[0 12],'Box','off','TickDir','out','FontSize',8,...
    'LineWidth',0.35,'XTick',1:4:9,'YTick',0:4:8);
set(h(2:3),'YTickLabel',[]);
ylabel(h(1),'Lick rate (Hz)');
xlabel(h(1),'Session');

%% correlation b/w consummatory lick rate & dopamine response to reward
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 2.5 3.5]);
hold on;
bar(0.5,mean(rcorr),1,'FaceColor',[0.8 0.8 0.8],'EdgeColor','k','LineWidth',0.35);
scatter(rand(nMouse,1),rcorr,2,'k','filled');
errorbar(0.5,mean(rcorr),std(rcorr)/sqrt(nMouse),'k','Linewidth',0.5);
ylim([-0.3 0.3]);
xlim([-0.5 1.5])
set(gca,'Box','off','TickDir','out','FontSize',8,'LineWidth',0.35,...
    'XTick',[],'YTick',-0.3:0.3:0.3)
ylabel('r (DA vs. consumption lick rate)');

%% doapmine response around consummatory and nonconsummatory lick onset
clr = {[0.6 0.6 0.6],[0 0 0];...
    [0.6 0.6 1],[0 0 1]};

fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 4 3.5]);
hold on;
for i = 1:2
    if i==1
        data = squeeze(nanmean(avepsth_nonconsum,2));
    else
        data = squeeze(nanmean(avepsth_consum,2));
    end
    plot(time/1000,data,'Color',clr{i,1},'LineWidth',0.35);
    plot(time/1000,nanmean(data),'Color',clr{i,2},'LineWidth',1);
end
xlim([-1.5 3]);
plot([0 0],[-2 6],'k:');
set(gca,'Box','off','XTick',-1.5:1.5:3,'YTick',-2:2:6,'TickDir','out',...
    'FontSize',8,'LineWidth',0.35,'YLim',[-2 6]);
xlabel('Time from lick bout onset (s)');
ylabel('Normalized DA response');
