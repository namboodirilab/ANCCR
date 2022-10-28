clearvars; clc; 

rng(7)

directory = 'D:\OneDrive - University of California, San Francisco\Huijeong\DA';
% directory = 'D:\OneDrive - UCSF\Photometry';

mouseList = {'HJ_FP_M2';'HJ_FP_M3';'HJ_FP_M4';'HJ_FP_F1';'HJ_FP_F2';'HJ_FP_M6';'HJ_FP_M7'};
nMouse = length(mouseList);

startday = [29,24,32,19,24,19,20; 70,61,66,42,46,47,43];
aucdata = cell(nMouse,2);
win = [3000 4000; 9000 10000];
for iM = 1:nMouse
    iM
    for i = 1:2
        behfile = FindFiles('Events_cues.mat','StartingDirectory',[directory,'\',mouseList{iM}],'CheckSubdir',1);
        out = ~cellfun(@(x) contains(x,'Day'),behfile);
        behfile(out) = [];
        temp = cellfun(@(x) strsplit(x,'Day'),behfile,'UniformOutput',false);
        temp = cellfun(@(x) strsplit(x{2},{'\','_'}),temp,'UniformOutput',false);
        day = cellfun(@(x) str2double(x{1}),temp);
        if i==1
            behfile = behfile(day==startday(i,iM));
        else
            [day,sortidx] = sort(day);
            behfile = behfile(sortidx);
            behfile = behfile(find(day<=startday(i,iM),2,'last'));
        end
        load(behfile{1});
        CSrw = unique(CS(fxreward==1));
        if i==2
            behfile(1) = [];
            load(behfile{1});
        end
        
        load([fileparts(behfile{1}),'\Photometry.mat']);
        
        [timecs,cssignal] = alignsignal2event(T(:,1),dff',CStime(CS==CSrw),[-2000 12000],10);
        aucbase = aucsignal(cssignal,timecs,[-1000 0]);
        auccs = aucsignal(cssignal,timecs,win(i,:))/diff(win(i,:));
        
        aucdata{iM,i} = auccs-aucbase;
    end
end

%%
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
%%
cd('D:\heejeong\OneDrive - University of California, San Francisco\figures\manuscript\dopamine_contingency\figS8');
print(fHandle,'-depsc','-painters','omission_rsp.ai');

