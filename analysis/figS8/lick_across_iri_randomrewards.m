clearvars; clc; close all;

directory = 'D:\OneDrive - UCSF\Huijeong\DA\DB_longITI_C1_';
mouseList = {'F1';'F2';'M1';'M2';'M3';'M4'};
nMouse = length(mouseList);

lickrate = nan(nMouse,6,2);
for iM = 1:nMouse
    behfile = findfiles('Events_randomrewards.mat',[directory,mouseList{iM}],1,'day');
    iri= cellfun(@(y) str2double(y(6)),cellfun(@(x) strsplit(fileparts(x),{'day','_'}),...
        behfile,'UniformOutput',false));
    [~,idx] = sort(iri);
    behfile = behfile(idx);
    nFile = length(behfile);
    
    for iF = 1:nFile
        load(behfile{iF});
        lickrate(iM,iF,2) = sum(~lickconsumidx)/((sessionendtime-nansum(consumboutlength))/10^3);
        lickrate(iM,iF,1) = sum(lickconsumidx)/(nansum(consumboutlength)/10^3);
    end
end


%%
avelickrate = movmean(lickrate,2,2);
avelickrate = avelickrate(:,2:2:end,:);
titlelist = {'consummatory';{'non';'comsummatory'}};
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 4 4]);
for i = 1:2
   h(i) = subplot(1,2,i); hold on;
   plot([6,12,30],avelickrate(:,:,i),'Color',[0.6,0.6,0.6])
   errorbar([6,12,30],mean(avelickrate(:,:,i),1),std(avelickrate(:,:,i),[],1)/sqrt(nMouse),'k','CapSize',3)
   title(titlelist{i})
end
set(h,'XTick',[6 12 30],'Box','off','TickDir','out',...
    'FontSize',8,'LineWidth',0.35,'XLim',[1 35]);
set(h(1),'YTick',4:7,'YLim',[4 7]);
set(h(2),'YTick',0:0.5:1.5,'YLim',[0 1.5]);
ylabel(h(1),'Lick rate (Hz)');
xlabel('IRI (s)')