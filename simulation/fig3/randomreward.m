clearvars; clc;
rng(7)

%% parameter set up
% task parameters
meanITI = 12;
numrewards = [50000, 20000, 2000];

% anccr model parameters
samplingperiod = 0.2;
alpha_anccr.exponent = 0.1;
alpha_anccr.init = 0.25;
alpha_anccr.min = 0.02;        
alpha_r = 0.2;
w = 0.5;               
k = 1;                 
minimumrate = 10^(-3);
maximumjitter = 0.1;
beta = [1];
threshold = 0.6;
T = meanITI*1.2;

% rpe model parameters - csc/microstimulus
alpha_rpe = [0.05, 0.02];
gamma = [0.95, 0.98];
lambda = [0, 0.95];
statesize = 0.2;
sigma = 0.08;
nmicrostimulus = 20;
d = 0.99;

nIter = 100;
rwrsp = cell(3,1);
iri = cell(3,1);
%% simulation
for iiter = 1:nIter
    iiter
    eventlog = simulateBackgroundRewards(numrewards(1),meanITI,1,1,0);
    eventlog(:,2) = eventlog(:,2)+1;
    for imdl = 1:3
        if iiter==1
            [rwrsp{imdl},iri{imdl}] = deal(nan(numrewards(imdl),nIter));
        end
        switch imdl
            case 1
                [DA,valuetimeline,eventtimeline,statetimeline] =...
                    simulateCSC(eventlog(1:numrewards(imdl),:),1,statesize,alpha_rpe(imdl),...
                    gamma(imdl),lambda(1),[],1,nan,meanITI*3);
            case 2
                [DA,~,eventtimeline] = simulatemicrostimulus(eventlog(1:numrewards(imdl),:),...
                    1,nmicrostimulus,statesize,alpha_rpe(imdl),gamma(imdl),lambda(2),sigma,d,[],meanITI*3);
            case 3
                DA = calculateANCCR(eventlog(1:numrewards(imdl),:),T,alpha_anccr,k,samplingperiod,w,threshold,...
                    minimumrate,beta,alpha_r,maximumjitter);
                eventtimeline = eventlog(1:numrewards(imdl),1);
        end
        rwrsp{imdl}(1:sum(eventtimeline(:,1)==1),iiter) = DA(eventtimeline(:,1)==1);
        if imdl<3
            iri{imdl}(1:sum(eventtimeline(:,1)==1),iiter) = [nan;diff(find(eventtimeline(:,1)==1))*statesize];
        else
            iri{imdl}(:,iiter) = [nan;diff(eventlog(1:numrewards(imdl),2))];
        end
    end
end
%%
dir = 'D:\OneDrive - University of California, San Francisco\figures\manuscript\dopamine_contingency\revision';
save([dir,'\data\randomrewards.mat'],'iri','rwrsp');
%%
dir = 'D:\OneDrive - University of California, San Francisco\figures\manuscript\dopamine_contingency\revision';
load([dir,'\data\randomrewards.mat'],'iri','rwrsp');
%%
ft = fittype('A*(1-2^(-(x/L)^S))+b'); %Weibull function; asymptote (b+A), latency (L), abruptness (S)

asymptotetrial = nan(3,1);
for imdl = 1:3
    numrwstates = find(sum(isnan(rwrsp{imdl}),2)==0,1,'last');
    [curve1,gof1] = fit([1:numrwstates]',nanmean(rwrsp{imdl}(1:numrwstates,:),2),ft,...
        'Lower',[-2 10 1 1],'Upper',[2 numrwstates 2 3]);
    [curve2,gof2] = fit([1:numrwstates]',nanmean(rwrsp{imdl}(1:numrwstates,:),2),ft,...
        'Lower',[-4 0 -1 -2],'Upper',[4 numrwstates 2 1]);
    if gof1.rsquare>gof2.rsquare
        curve = curve1;
    else
        curve = curve2;
    end
    coef = coeffvalues(curve);
    asymptote = coef(1)+coef(4);
    if coef(4)>coef(1)
        asymptotetrial(imdl) = find(curve(1:numrwstates)>coef(1)*0.95+coef(4),1,'last');
    else
        asymptotetrial(imdl) = find(curve(1:numrwstates)<coef(1)*0.95+coef(4),1,'last');
    end
    subplot(1,3,imdl)
    hold on;
    plot([1:numrwstates]',nanmean(rwrsp{imdl}(1:numrwstates,:),2));
    plot(curve)
end



for imdl = 1:3
    r(:,imdl) = cellfun(@(x) corr([1:asymptotetrial(imdl)]',x),...
        mat2cell(rwrsp{imdl}(1:asymptotetrial(imdl),:),asymptotetrial(imdl),ones(1,nIter)));
end

%%
model = {'RPE (CSC)';'RPE (MS)';'ANCCR'};
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 3 4]);
hold on;
for imdl = 1:3
    if imdl<3
       x = 1.5*imdl-1;
    else
        x = 1.5*imdl;
    end
    r = cellfun(@(x) corr([1:asymptotetrial(imdl)]',x),...
        mat2cell(rwrsp{imdl}(1:asymptotetrial(imdl),:),asymptotetrial(imdl),ones(1,nIter)));
   bar(x,mean(r),1,'FaceColor',[0.8 0.8 0.8],'EdgeColor','k','LineWidth',0.35);
   errorbar(x,mean(r),std(r)/sqrt(nIter),'k','Linewidth',0.5,'CapSize',3);
   [~,p,~,stat] = ttest(r)

end
ylim([-0.3 0.5])
set(gca,'Box','off','TickDir','out','FontSize',7,'LineWidth',0.35,...
    'XTick',[0.5 2 4.5],'XTickLabel',model,'XTickLabelRotation',45,...
    'YTick',-0.4:0.2:0.4,'XLim',[-0.5 5.5]);
ylabel('r (predicted DA vs. trial)','FontSize',8);
print(fHandle,'-depsc','-painters',[dir,'\randomrewards_trial_bar.ai']);

%%
imdl = 2;
iiter = 7;
% in = 1:size(DArsp,1)<lasttrial(imdl) & iri(:,iiter)'>3;
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 3.1 3.3]);
hold on;
[r,p] = corr([1:asymptotetrial(imdl)]',rwrsp{imdl}(1:asymptotetrial(imdl),iiter));
scatter(find(1:asymptotetrial(imdl)),rwrsp{imdl}(1:asymptotetrial(imdl),iiter),2,'k','filled');
beta = glmfit(find(1:asymptotetrial(imdl)),rwrsp{imdl}(1:asymptotetrial(imdl),iiter));
plot([0 asymptotetrial(imdl)],beta(1)+beta(2)*[0 asymptotetrial(imdl)],'r','LineWidth',1);
text(asymptotetrial(imdl)*0.1,2,['r = ',num2str(round(r*100)/100)],'FontSize',7);
text(asymptotetrial(imdl)*0.1,1.8,['p = ',num2str(round(p*100)/100)],'FontSize',7);
xlim([0 asymptotetrial(imdl)]);
set(gca,'Box','off','TickDir','out','FontSize',8,'LineWidth',0.35,...
    'YTick',0.5:0.5:2,'XTick',0:500:asymptotetrial(imdl),'YTick',0:1:3);
xlabel('Trial','FontSize',8);
ylabel('DA response','FontSize',8);
title('Microstimulus');
% ylim([0.25 2.5]);
print(fHandle,'-depsc','-painters',[dir,'\randomrewards_trial_microstimulus.ai']);



%% 
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 3 4]);
hold on;
for imdl = 1:3
    if imdl<3
       x = 1.5*imdl-1;
    else
        x = 1.5*imdl;
    end
    r = cellfun(@(x,y,z) corr(y(x>3),x(x>3),'rows','complete'),...
        mat2cell(iri{imdl},numrewards(imdl),ones(nIter,1)),...
        mat2cell(rwrsp{imdl},numrewards(imdl),ones(nIter,1)));
   bar(x,mean(r),1,'FaceColor',[0.8 0.8 0.8],'EdgeColor','k','LineWidth',0.35);
   errorbar(x,mean(r),std(r)/sqrt(nIter),'k','Linewidth',0.5,'CapSize',3);
      [~,p,~,stat] = ttest(r)
end
ylim([-0.8 0.8])
set(gca,'Box','off','TickDir','out','FontSize',7,'LineWidth',0.35,...
    'XTick',[0.5 2 4.5],'XTickLabel',model,'XTickLabelRotation',45,...
    'YTick',-0.8:0.4:0.8,'XLim',[-0.5 5.5]);
ylabel('r (predicted DA vs. IRI)','FontSize',8);
cd(dir);
print(fHandle,'-depsc','-painters',[dir,'\randomrewards_iri_bar.ai']);

%%
iiter = 2;
imdl = 2;
in = iri{imdl}(:,iiter)>3;
[r,p] = corr(iri{imdl}(in,iiter),rwrsp{imdl}(in,iiter));

fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 3.1 3.3]);
hold on;
scatter(iri{imdl}(in,iiter),rwrsp{imdl}(in,iiter),2,'k','filled');
beta = glmfit(iri{imdl}(in,iiter),rwrsp{imdl}(in,iiter));
plot([0 max(iri{imdl}(in,iiter))],beta(1)+beta(2)*[0 max(iri{imdl}(in,iiter))],'r');
text(50,2.3,['r = ',num2str(round(r*100)/100)],'FontSize',7);
text(50,2,['p = ',num2str(p)],'FontSize',7);
xlim([0 max(iri{imdl}(in,iiter))]);
ylim([0 2.5])
set(gca,'Box','off','TickDir','out','FontSize',8,'LineWidth',0.35,...
    'XTick',0:50:max(iri{imdl}(in,iiter)),'YTick',0:2);
xlabel('IRI','FontSize',8);
ylabel('Predcted DA response','FontSize',8);
title('Microstimulus');
print(fHandle,'-depsc','-painters',[dir,'\randomrewards_iri_microstimulus.ai']);





