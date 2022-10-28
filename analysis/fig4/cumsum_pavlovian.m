clearvars; clc; close all;

rng(7)

directory = 'D:\OneDrive - University of California, San Francisco\Huijeong\DA';

mouseList = {'HJ_FP_M2';'HJ_FP_M3';'HJ_FP_M4';'HJ_FP_F1';'HJ_FP_F2';'HJ_FP_M6';'HJ_FP_M7'};
nMouse = length(mouseList);

% initial learning, cue duration change, extinction
startday = [16,29,70,38; 12,24,61,28; 13,32,66,38;...
    7,19,42,23; 13,24,46,28; 7,19,47,23; 6,20,43,24]; % first day or each condition
endday = [28,32,72,43; 23,27,63,38; 31,37,67,47;...
    18,22,43,27; 23,27,47,33; 18,22,50,30; 19,23,44,29]; % last day of each condition
learnedday = [22,17,25,12,17,12,14]; % first day w/ anticipatory behavior

ndays = [endday-startday]+1;
ndaybeforestart = [0,3,1,1];
ratio = [1.5,1.9,2.5,2];

[normcum_lick_cs,normcum_lick_base,normcum_lick,normcum_auc] = deal(cell(nMouse,4));
[abruptness_lick,abruptness_auc,changetrial_lick,changetrial_auc] = deal(nan(nMouse,4));
[daic_lick,daic_auc] = deal(nan(nMouse,2));

for iM = 1:nMouse
    behfile = findfiles('Events_cues.mat',[directory,'\',mouseList{iM},'\Pavlovian'],1,'Day');
    days = cellfun(@(y) str2double(y(4)),cellfun(@(x) strsplit(fileparts(x),{'Day','_'}),...
        behfile,'UniformOutput',false));
    [days,sortidx] = sort(days);
    behfile = behfile(sortidx);
    
    for iC = 1:4
        %% load data for each condition &  calculate auc, lick number
        in = find(days>=startday(iM,iC)-ndaybeforestart(iC) & days<=endday(iM,iC));
        nin = length(in);
        [lickcs,lickbase,auc] = deal(cell(nin,1));
        for i = 1:nin
            load(behfile{in(i)});
            try
                load([fileparts(behfile{in(i)}),'\Photometry.mat']);
            catch
                continue;
            end
            
            CSrw = unique(CS(fxreward==1));
            
            [timecs,cssignal] = alignsignal2event(T(:,1),dff',CStime(CS==CSrw),[-2000 10000],10);
            aucbase = aucsignal(cssignal,timecs,[-1000 0]);
            auccs = aucsignal(cssignal,timecs,[0 2000])/2;
            auc{i} = auccs-aucbase;
            
            if iC==2
                % median ancitipatory lick time
                lickcs{i} = cellfun(@(x) nanmedian(licktime(licktime>=x & licktime<=x+csduration+1000))-x,num2cell(CStime(CS==CSrw)));
                out = isnan(lickcs{i}); % exclude trial w/o anticipatory lick
                lickcs{i} = lickcs{i}(~out);
                auc{i} = auc{i}(~out);
            else
                lickcs{i} = numevent(licktime,CStime(CS==CSrw),[0 csduration+1000],1);
                lickbase{i} = numevent(licktime,CStime(CS==CSrw),[-1000 0],1);
            end
        end
        
        if iC<4
            %% make cumsum plot, calculate abruptness and change trial from that
            % according to number of trials before condition change, determine
            % total number of trials to use; this is to align data across animals
            if iC==1
                inbefore = days(in)<learnedday(iM);
            else
                inbefore = days(in)<startday(iM,iC);
            end
            ntrialbefore = length(cell2mat(lickcs(inbefore)));
            ntotaltrial = round(ntrialbefore*ratio(iC));
            
            lickcum_cs = cumsum(cell2mat(lickcs));
            if iC==2
                lickcum = lickcum_cs;
            else
                lickcum_base = cumsum(cell2mat(lickbase));
                normcum_lick_cs{iM,iC} = lickcum_cs(1:ntotaltrial)/lickcum_cs(ntotaltrial);
                normcum_lick_base{iM,iC} = lickcum_base(1:ntotaltrial)/lickcum_cs(ntotaltrial);
                lickcum = lickcum_cs-lickcum_base;
            end
            auccum = cumsum(cell2mat(auc));
            
            normcum_lick{iM,iC} = lickcum(1:ntotaltrial)/lickcum(ntotaltrial);
            normcum_auc{iM,iC} = auccum(1:ntotaltrial)/auccum(ntotaltrial);
            
            lickcum_s = conv(lickcum(1:ntotaltrial),fspecial('Gaussian',[1 5*2],2),'valid'); % smooth graph
            auccum_s = conv(auccum(1:ntotaltrial),fspecial('Gaussian',[1 5*2],2),'valid');
            lickcum_s = lickcum_s/lickcum_s(end); % normalized cumsum plots
            auccum_s = auccum_s/auccum_s(end);
            
            % distance from diagonal
            d_lick = point2line([[5:ntotaltrial-5]',lickcum_s],[5,lickcum_s(1)], [ntotaltrial-5,1]);
            d_auc = point2line([[5:ntotaltrial-5]',auccum_s],[5,auccum_s(1)], [ntotaltrial-5,1]);
            
            % find abruptness & change trial
            [abruptness_lick(iM,iC),idx_lick] = max(d_lick);
            [abruptness_auc(iM,iC),idx_auc] = max(d_auc);
            changetrial_lick(iM,iC) = idx_lick;
            changetrial_auc(iM,iC) = idx_auc;
            
            
            %% testing changing vs. constant model for cue duration change experiment
            if iC==2
                lickdata = cell2mat(lickcs);
                lickdata = lickdata(1:ntotaltrial);
                aucdata = cell2mat(auc);
                aucdata = aucdata(1:ntotaltrial);
                
                for imdl = 1:2
                    if imdl==1
                        ft = fittype('A*(1-2^(-(x/L)^S))+b'); %Weibull function; asymptote (A+b), latency (L), abruptness (S)
                        [curve_lick,gof_lick] = fit([1:ntotaltrial]',lickdata,ft,...
                            'Lower',[0 ntrialbefore-10 0 0],'Upper',[9000 ntotaltrial 100 4000]);
                        [curve_auc,gof_auc] = fit([1:ntotaltrial]',aucdata,ft,...
                            'Lower',[-2 ntrialbefore-10 0 -10],'Upper',[50 ntotaltrial 10 10]);
                    else
                        ft_constant = fittype('0*x+b');
                        [curve_lick,gof_lick] = fit([1:ntotaltrial]',lickdata,ft_constant);
                        [curve_auc,gof_auc] = fit([1:ntotaltrial]',aucdata,ft_constant);
                    end
                    rss_lick = sum((lickdata-curve_lick(1:ntotaltrial)).^2);
                    rss_auc = sum((aucdata-curve_auc(1:ntotaltrial)).^2);
                    daic_lick(iM,imdl) = 2*length(coeffvalues(curve_lick))+ntotaltrial*log(rss_lick);
                    daic_auc(iM,imdl) = 2*length(coeffvalues(curve_auc))+ntotaltrial*log(rss_auc);
                end
            end
        end
    end
end


for iC = 1:4
    if iC<4
        %% plotting normalized cumsum plot
        fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 5.5 3.3]);
        for i = 1:2
            axes('Position',axpt(2,1,i,1,axpt(10,10,1:10,2:8),[0.15 0.05]));
            hold on;
            if i==1
                data = normcum_lick(:,iC);
            else
                data = normcum_auc(:,iC);
            end
            for iM = 1:nMouse
                plot([1:length(data{iM})]/length(data{iM}),data{iM},'Color',[0.6,0.6,0.6],'LineWidth',0.35);
            end
            plot([0,1],[0,1],'k--');
            set(gca,'YLim',[-0.2 1.05],'Box','off','TickDir','out','FontSize',7,'LineWidth',0.35,'XTick',0:0.5:1,'YTick',0:0.5:1);
            if i==1
                ylabel({'Normalized';'cumsum (lick rate)'});
            else
                ylabel({'Normalized cumsum';'(DA cue response)'});
            end
            if iC>1
                plot(repmat(1/ratio(iC),1,2),[-0.2 1.05],'k:','LineWidth',0.35);
            end
        end
        
        %% change trial / abruptness / aic
        fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 4.5 3.3]);
        for i = 1:2
            subplot(1,2,i);
            hold on;
            if i==1
                data = [abruptness_lick(:,iC),abruptness_auc(:,iC)];
            else
                if iC==2
                    data = [diff(daic_lick,[],2),diff(daic_auc,[],2)];
                    plot([0.5 2.5],[0 0],'k:','LineWidth',0.35);
                else
                    data = [changetrial_lick(:,iC),changetrial_auc(:,iC)];
                end
            end
            plot([1,2],data,'Color',[0.6 0.6 0.6],'LineWidth',0.35);
            errorbar(1:2,mean(data),std(data)/sqrt(nMouse),'k');
            set(gca,'Box','off','TickDir','out','FontSize',8,'LineWidth',0.35,...
                'XTick',1:2,'XTickLabel',{'Lick','DA'},'XTickLabelRotation',45,'XLim',[0.5 2.5]);
            if i==1
                set(gca,'YTick',0:0.5:1,'YLim',[0,1]);
                ylabel('Abruptness of change');
            else
                if iC==2
                    ylabel('\DeltaAIC(constant-changing)');
                else
                    ylabel('Change trial');
                end
            end
        end
    end
    %% plotting unnormalized cumsum plot (figS9)
    clr = {[0.6 0.6 0.6],[0 0 0],[1 0 0]};
    fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 18 3.3]);
    for iM = 1:nMouse
        axes('Position',axpt(nMouse,1,iM,1,axpt(1,10,1,1:8),[0.02 0.05]));
        hold on;
        for i = 1:3
            if iC==2 & i==1
                continue;
            end
            switch i
                case 1
                    data = normcum_lick_base{iM,iC};
                case 2
                    if iC==2
                        data = normcum_lick{iM,iC};
                    else
                        data = normcum_lick_cs{iM,iC};
                    end
                case 3
                    data = normcum_auc{iM,iC};
            end
            plot(data,'Color',clr{i},'LineWidth',1);
            
        end
        plot([0,length(data)],[0,1],'k--');
        set(gca,'YLim',[-0.2 1.05],'Box','off','TickDir','out','FontSize',7,'LineWidth',0.35,'YTick',0:0.5:1);
        if iM==1
            ylabel({'Normalized';'cumsum'});
        else
            set(gca,'YTickLabel',[]);
        end
        xlabel('Trial');
        if iC>1
            plot(repmat(length(data)/ratio(iC),1,2),[-0.2 1.05],'k:','LineWidth',0.35);
        end
    end
end
