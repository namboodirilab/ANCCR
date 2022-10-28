clearvars; close all; clc;
rng(2)

% task parameter
numcue = 800;
cuerewdelay = 3;
cuecuedelay = 1.5;
consumdelay = 3;
meanITI = 60;
inhibitionstrength = -0.6;

% model parameter
alpha = 0.05; % learning rate
gamma = 0.95; % discount factor
lambda = 0;
numstimulus = 20; % number of microstimulus
d = 0.99; % memory decay constant
statesize = 0.2;
sigma = 0.08; % width of basis function

%%

nIter = 20;
darsp = cell(3,3); % row: different inh levels; column: different events
for iIter = 1:nIter
    eventlog = simulateEventChain(numcue, 2, 4, meanITI, ...
        meanITI*3, cuecuedelay, cuerewdelay, 1, consumdelay);
    for iinh = 1:3 % constant inh; no inh; increasing inh
        printf('%d iteration: %d inh level\n',iIter,iinh);
        switch iinh
            case 1
                inhibitionlog = [eventlog(eventlog(:,1)==2,2),...
                eventlog(eventlog(:,1)==2,2)+cuecuedelay,...
                ones(numcue,1)*inhibitionstrength];
            case 2
                inhibitionlog = [];
            case 3
                inhibitionlog = [eventlog(eventlog(:,1)==2,2),...
                eventlog(eventlog(:,1)==2,2)+cuecuedelay,...
                ones(numcue,1)*-0.1, ones(numcue,1)*inhibitionstrength];
        end
                [rpetimeline,valuetimeline,eventtimeline,~,inhibitiontimeline] =...
                    simulateCSC(eventlog,3,statesize,alpha,gamma,lambda,inhibitionlog);
            
        for ie = 1:3 % 1:cs1, 2:cs2, 3:reward
            inevent = eventtimeline(:,1)==ie;
            darsp{iinh,ie}(:,iIter) = rpetimeline(inevent);
            if ie==1
                valuemap{iinh}(:,:,iIter) =...
                    cell2mat(cellfun(@(x) valuetimeline([-1/statesize:5/statesize]+x)',...
                    num2cell(find(inevent)),'UniformOutput',false));
            end
        end
    end
end
%%
close all
clr_light = [0.6,0.6,0.6; 0.6,0.6,1];
cuelist = {'CS1,','CS2'};
clr = [0 0 0; 0 0 1];

dir = 'D:\OneDrive - UCSF\figures\manuscript\dopamine_contingency\revision\fig6_new';
% dir = 'D:\OneDrive - University of California, San Francisco\figures\manuscript\dopamine_contingency\revision\fig6_new';
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 4 3]);
axes('Position',axpt(5,5,2:5,1:4))
hold on;
for i= 1:2
% s = std(darsp{2,i},[],2)'/sqrt(nIter);
plot(1:numcue,darsp{2,i},'Color',clr_light(i,:),'LineWidth',0.35)
% fill([1:numcue, flip(1:numcue)],[m+s, flip(m-s)],clr_light(i,:),'EdgeColor','none')
% plot(1:numcue,m,'Color',clr(i,:),'LineWidth',1)
plot([20 60],[1.1 1.1]-i*0.1,'Color',clr(i,:),'LineWidth',1);
text(80,1.1-i*0.1,cuelist{i})
end
for i = 1:2
    m = mean(darsp{2,i},2)';
    plot(1:numcue,m,'Color',clr(i,:),'LineWidth',1)

end
xlabel('Trial')
ylabel('Predicted DA response')
set(gca,'FontSize',8,'Box','off','TickDir','out','LineWidth',0.35,'XTick',0:200:800,...
    'YTick',-0.1:0.1:0.4,'YLim',[-0.05 0.4],'XTickLabelRotation',45)
print(fHandle,'-depsc','-painters',[dir,'\sequential_csc.ai'])

%% 
darsp2 = darsp;
clr ={'r','k'};
yaxisbreak = [-5.1 -0.1];

fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 3 3]);
axes('Position',axpt(5,5,3:5,1:4))
hold on;
for imdl = 1:2
    if imdl==2
        load(['sequential_anccr.mat']);
    else
        load(['sequential_rpe.mat']);
    end
    data = cell2mat(cellfun(@(x) mean(x(end-100:end,1:nIter),1),darsp(:,1),'UniformOutput',false))';
    data(data<yaxisbreak(1)) = data(data<yaxisbreak(1))+diff(yaxisbreak);
    for i = 1:2
        scatter(rand(nIter,1)+1.5*(imdl-1),data(:,i),2,clr{i},'filled');
        errorbar(0.5+1.5*(imdl-1),mean(data(:,i)),std(data(:,i))/sqrt(nIter),'LineWidth',0.5,'Color',clr{i});
    end
    [~,p,~,stat] = ttest2(data(:,1),data(:,2))
end
plot([-0.5, 3],[0 0],'k:','LineWidth',0.35);
xlim([-0.5 3]);
ylim([-2 1.5]);
yticks = -2:1;
yticklabels = yticks;
yticklabels(yticklabels<yaxisbreak(2)) = yticklabels(yticklabels<yaxisbreak(2))-diff(yaxisbreak);
set(gca,'Box','off','TickDir','out','FontSize',8,'LineWidth',0.35,'XTick',[0.5,2],...
    'XTickLabel',{'Model 1','Model 2'},'YTick',yticks,'YTickLabel',yticklabels,'XTickLabelRotation',45);
ylabel('Predicted DA response');
print(fHandle,'-depsc','-painters',[dir,'\sequential_csc_cs1.ai'])

%%
close all
darsp = darsp2;
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 3 3.5]);
axes('Position',axpt(5,5,3:5,1:4))
hold on;
for imdl = 1:2
    if imdl==2
        load([dir,'\sequential_anccr.mat']);
    end
    data = cell2mat(cellfun(@(x) mean(x(end-50:end,1:nIter),1),darsp(2,1:2),'UniformOutput',false)')';
    data = data(:,2)./data(:,1);
    bar(0.5+1.5*(imdl-1),mean(data),1,'FaceColor',[0.6 0.6 0.6])
    scatter(rand(nIter,1)+1.5*(imdl-1),data,2,'k','filled')
    errorbar(0.5+(imdl-1)*1.5,mean(data),std(data)/sqrt(nIter),'k','LineWidth',0.5)
end
set(gca,'Box','off','TickDir','out','FontSize',8,'LineWidth',0.35,'XTick',[0.5,2],'XTickLabel',...
    {'Model1','Model2'},'YTick',0:0.5:1,'XTickLabelRotation',45)
ylim([-0.1 1])
ylabel({'Norm. CS2 response';'(1 = CS1 response)'})
print(fHandle,'-depsc','-painters',[dir,'\sequential_csc_cs2_last50.ai'])

