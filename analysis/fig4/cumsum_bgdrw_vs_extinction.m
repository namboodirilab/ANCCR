clearvars; clc; close all;

rng(7)

directory = 'D:\OneDrive - UCSF\Huijeong\DA';

mouseList = {'HJ_FP_M2';'HJ_FP_M3';'HJ_FP_M4';'HJ_FP_F1';'HJ_FP_F2';'HJ_FP_M6';'HJ_FP_M7'};
nMouse = length(mouseList);



% background reward sessions, extinction sessions
startday = [33,70; 28,61; 38,66; 23,42; 28,46; 23,47; 24,43]; % first day or each condition
endday = [43,72; 38,63; 47,67; 27,43; 33,47; 30,50; 29,44]; % last day of each condition

ndaybeforestart = 2; % number of pre-condition sessions (for pooling only)
ntrialbeforestart = 50; % number of trials before start of each condition 
ntotaltrial = 140; % total number of trials 

normcum_auc = nan(nMouse,2,ntotaltrial);

for iM = 1:nMouse
    iM
    % pool sessions for each animal
    behfile = findfiles(mouseList{iM},[directory,'\',mouseList{iM},'\Pavlovian'],1,'Day','.mat');
    days = cellfun(@(y) str2double(y(4)),cellfun(@(x) strsplit(fileparts(x),{'Day','_'}),...
        behfile,'UniformOutput',false));
    [days,sortidx] = sort(days);
    behfile = behfile(sortidx);
    
    for iC = 1:2
        % find sessions in each condition
        in = find(days>=startday(iM,iC)-ndaybeforestart & days<=endday(iM,iC));
        nin = length(in);
        [lickcs,lickbase,auc] = deal(cell(nin,1));
        for i = 1:nin
            load(behfile{in(i)});
            try
                load([fileparts(behfile{in(i)}),'\Photometry.mat']);
            catch
                continue;
            end

            % load task related information
            eventtime = eventlog(:,2);
            eventindex = eventlog(:,1);
            nosolenoidflag = eventlog(:,3);
            [~,~,CS,CStime,~,~,fxreward] = eventfrompavlovian(eventtime,eventindex,nosolenoidflag);

            if iC==1
            CSrw = unique(CS(fxreward==1)); % find CS+ identity
            end
            
            % align dff to CS+ onset  
            [timecs,cssignal] = alignsignal2event(T(:,1),dff(:),CStime(CS==CSrw),[-2000 10000],10);
            
            % calculate AUC
            aucbase = aucsignal(cssignal,timecs,[-1000 0]); % AUC during pre-cue baseline
            auccs = aucsignal(cssignal,timecs,[0 2000])/2; % AUC during CS+
            
            % normalize CS+ response by baseline response
            auc{i} = auccs-aucbase; 
        end
        
        % trim data for population average later
        n = length(cell2mat(auc(1:ndaybeforestart)));
        data = cell2mat(auc);
        data = data([1:ntotaltrial]+n-ntrialbeforestart);

        % calculate normalized cumulative sum of dopamine response over 
        % trials across sessions
        normcum_auc(iM,iC,:) = cumsum(data)/sum(data(1:ntrialbeforestart));
    end
end

%%
% plot fig4P
clr = {[1 0.6 0.6],[1 0 0];[0.6 0.6 0.6],[0 0 0]};
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 3.5 3.2]);
hold on;
for iC = 1:2
    plot([1:ntotaltrial]/ntotaltrial,squeeze(normcum_auc(:,iC,:)),'Color',clr{iC,1},'LineWidth',0.35)
    plot([1:ntotaltrial]/ntotaltrial,mean(squeeze(normcum_auc(:,iC,:)),1),'Color',clr{iC,2},'LineWidth',1)
end
plot(repmat(ntrialbeforestart/ntotaltrial,1,2),[0 2.5],'k:','LineWidth',0.35);
xlim([0 1]);
ylim([0 2.5]);
xlabel('Normalized trial');
ylabel({'Normalized';'cumsum (DA response)'});
set(gca,'Box','off','TickDir','out','FontSize',8,'LineWidth',0.35,...
    'XTick',[0 0.5 1],'YTick',0:2);