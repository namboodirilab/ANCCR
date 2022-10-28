clearvars; clc; close all;

rng(7)

% directory = 'D:\OneDrive - UCSF\Huijeong\DA';
directory = 'D:\OneDrive - University of California, San Francisco\Huijeong\DA';
mouseList = {'HJ_FP_M2';'HJ_FP_M3';'HJ_FP_M4';'HJ_FP_F1';'HJ_FP_F2';'HJ_FP_M6';'HJ_FP_M7'};
nMouse = length(mouseList);

startday = [77,67,70,46,51,53,48];
endday = [78,68,72,49,54,55,50];
ndays = endday-startday+1;

[avecuersp,averwrsp] = deal(nan(max(ndays)+1,nMouse,2));
for iM = 1:nMouse
    iM
    behfile = findfiles('Events_poisson.mat',[directory,'\',mouseList{iM},'\Poisson'],1,'Day');
    days = cellfun(@(y) str2double(y(4)),cellfun(@(x) strsplit(fileparts(x),{'Day','_'}),...
        behfile,'UniformOutput',false));
    [days,sortidx] = sort(days);
    behfile = behfile(sortidx);
    in = find(days>=startday(iM) & days<=endday(iM));
    days = days(in);
    behfile = behfile(in);
    nin = length(in);

    for i = 1:nin
        load(behfile{i});
        try
            load([fileparts(behfile{i}),'\Photometry.mat']);
        catch
            continue;
        end
        CSrw = unique(CS(fxreward==1));
        cuerewdelay = mean(fxrwtime-CStime);

        [timecs,cssignal] = alignsignal2event(T(:,1),dff(:),CStime(CS==CSrw),[-2000 10000],0);
        aucbase = aucsignal(cssignal,timecs,[-500 0]);
        auccs = aucsignal(cssignal,timecs,[0 500]);
        
        [timerw,rwsignal] = alignsignal2event(T(:,1),dff(:),firstfxlicktime(CS==CSrw),[-2000 5000],10);
        aucrw = aucsignal(rwsignal,timerw,[0 500],isnan(firstfxlicktime(CS==CSrw)));

        % finding intermediate cues and previous cues of them
        intermediateidx = [false;diff(CStime)<cuerewdelay];
        trialidx = [[intermediateidx(2:end);0],intermediateidx];

        % if the delays from previous cue and to next cue are both shorter
        % than cuerewdelay, excluded all three cues from analysis
        trialidx(find(sum(trialidx,2)==2)-1,:) = 0;
        trialidx(find(sum(trialidx,2)==2)+1,:) = 0;
        trialidx(sum(trialidx,2)==2,:) = 0;

        % if delay b/w pair is shorter than 0.5s or intermediate cue to 
        % next reward delay is shorter than 0.5s excluded it from analysis
        shortintervalidx = find([false;diff(CStime)<=500] | [false;diff(CStime)>2500 & diff(CStime)<3000]);
        trialidx([shortintervalidx;shortintervalidx-1],:) = 0;
        trialidx = logical(trialidx);

        % for reward response, exclude rewards if animal didn't consumed it
        % before the next reward delivery
        inlick = [firstfxlicktime(1:end-1)<fxrwtime(2:end);0];

        aucbase(trialidx(:,2)) = aucbase(trialidx(:,1));
        auccs = auccs-aucbase;
        aucrw = aucrw-aucbase;

        for ii = 1:2
        avecuersp(i,iM,ii) = mean(auccs(trialidx(:,ii)));
        averwrsp(i,iM,ii) = mean(aucrw(trialidx(:,ii)&inlick));
        end
    end  
end
avecuersp = squeeze(nanmean(avecuersp,1));
averwrsp = squeeze(nanmean(averwrsp,1));

%%

data = avecuersp(:,2)./avecuersp(:,1);
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 2.2 2.9]);
hold on;
bar(0,mean(data),'FaceColor',[0.8 0.8 0.8],'EdgeColor','k','LineWidth',0.35);
errorbar(0,mean(data),std(data)/sqrt(nMouse),'k','Linewidth',0.5);
scatter(rand(nMouse,1)*0.8-0.4,data,2.5,'k','filled');
ylim([0 1]);
xlim([-1 1]);
set(gca,'Box','off','TickDir','out','FontSize',7,'LineWidth',0.35,'XTick',[],'YTick',0:0.5:1);
ylabel({'Ratio of cue response'; '(intermediate / previous)'},'FontSize',8)


%%
data = averwrsp(:,2)./averwrsp(:,1);
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 2.2 2.9]);
hold on;
bar(0,mean(data),'FaceColor',[0.8 0.8 0.8],'EdgeColor','k','LineWidth',0.35);
errorbar(0,mean(data),std(data)/sqrt(nMouse),'k','Linewidth',0.5);
scatter(rand(nMouse,1)*0.8-0.4,data,2.5,'k','filled');
plot([-1 1],[1 1],'k:','LineWidth',0.35);
ylim([0 3]);
xlim([-1 1]);
set(gca,'Box','off','TickDir','out','FontSize',7,'LineWidth',0.35,'XTick',[],'YTick',0:1:3);
ylabel({'Ratio of reward response'; '(intermediate / previous)'},'FontSize',8)
[~,p,~,stat] = ttest(data,1);
