%% 
clear all 
close all

% data
cd('\Users\rjc5t\Documents\Neuromechanics\DATA\Reward');
dat = readtable('Results\reward_table_long.csv');
dat = dat(dat.PREF == 0 & dat.TARG ~= 99 , :); % & dat.ID ~= 7
% aggregate data into means
aggdat = grpstats(dat, {'MASS','REWD'}, {'mean', 'sem'});

cd("C:\Users\rjc5t\Documents\Neuromechanics\DATA\UtilityModel\rtmtutility\")
load('simresults/simresults_young7.mat'); 
young = mysols;
% old
load('simresults/simresults_old7.mat'); 
old = mysols;


%% Define some variables
%range of alpha values - used in optimization
myalphas=20:100; % iterating over different levels of reward to see what happens

%range of SCALING on effort, alpha, and probability - used in optimization
myeffscales = 0.8:0.1:1.2;
myalphascales= 0.8:0.1:1.2;
myprobscales = 0.8:0.1:1.2;

% use this index to focus on alpha of a certain value
alphaind=find(myalphas==70);
alphascaleind=find(myalphascales==1);
effscaleind=find(myeffscales==1);
probscaleind= find(myprobscales==1);
myrange=find(myalphascales>=0.7);

% For plotting things
rwdlow = 35;
rwdhigh = 40;
alphalow = find(myalphas==rwdlow);
alphahigh = find(myalphas==rwdhigh);

% global parameters needed
global param
param.myc0to = -1; % Used in RT function
param.myc1to = 10; % Used in RT function
param.myc0 = -9; % MT young accuracy parameters; shifts logistic to the right with scaling
param.myc1 = 15; % MT young accuracy parameters; shifts logistic to the left with scaling


%% Process data

% Organize model data

% effort scale differences
for j=1:length(myeffscales)
    % young data
    effcurvesyoung=squeeze(young(:,j,alphascaleind,probscaleind,:));
    effcurvesyoung_mt(:,j)=effcurvesyoung(:,1);
    effcurvesyoung_rt(:,j)=effcurvesyoung(:,2);
    effsavings = ((effcurvesyoung_rt(alphahigh,j) - effcurvesyoung_rt(alphalow,j))+(effcurvesyoung_mt(alphahigh,j) - effcurvesyoung_mt(alphalow,j))); % total time savings
    effdeltasyoung(1,j) = (effcurvesyoung_mt(alphahigh,j) - effcurvesyoung_mt(alphalow,j));
    effdeltasyoung(2,j) = (effcurvesyoung_rt(alphahigh,j) - effcurvesyoung_rt(alphalow,j));
    effpropsyoung(1,j) = (effcurvesyoung_mt(alphahigh,j) - effcurvesyoung_mt(alphalow,j))/effsavings; % MT proportion of time savings
    effpropsyoung(2,j) = (effcurvesyoung_rt(alphahigh,j) - effcurvesyoung_rt(alphalow,j))/effsavings; % RT proportion of time savings
    
    % old data
    effcurvesold=squeeze(old(:,j,alphascaleind,probscaleind,:));
    effcurvesold_mt(:,j)=effcurvesold(:,1);
    effcurvesold_rt(:,j)=effcurvesold(:,2);
    effsavingsold = ((effcurvesold_rt(alphahigh,j) - effcurvesold_rt(alphalow,j))+(effcurvesold_mt(alphahigh,j) - effcurvesold_mt(alphalow,j))); % total time savings
    effdeltasold(1,j) = (effcurvesold_mt(alphahigh,j) - effcurvesold_mt(alphalow,j));
    effdeltasold(2,j) = (effcurvesold_rt(alphahigh,j) - effcurvesold_rt(alphalow,j));
    effpropsold(1,j) = (effcurvesold_mt(alphahigh,j) - effcurvesold_mt(alphalow,j))/effsavingsold; % MT proportion of time savings
    effpropsold(2,j) = (effcurvesold_rt(alphahigh,j) - effcurvesold_rt(alphalow,j))/effsavingsold; % RT proportion of time savings
end

% reward scale differences
for k=1:length(myalphascales)
    alphacurvesyoung=squeeze(young(:,effscaleind,k,probscaleind,:));
    alphacurvesyoung_mt(:,k)=alphacurvesyoung(:,1);
    alphacurvesyoung_rt(:,k)=alphacurvesyoung(:,2);
    rwdsavings = ((alphacurvesyoung_rt(alphahigh,k) - alphacurvesyoung_rt(alphalow,k))+(alphacurvesyoung_mt(alphahigh,k) - alphacurvesyoung_mt(alphalow,k)));
    rwddeltasyoung(1,k) = (alphacurvesyoung_mt(alphahigh,k) - alphacurvesyoung_mt(alphalow,k));
    rwddeltasyoung(2,k) = (alphacurvesyoung_rt(alphahigh,k) - alphacurvesyoung_rt(alphalow,k));
    rwdpropsyoung(1,k) = (alphacurvesyoung_mt(alphahigh,k) - alphacurvesyoung_mt(alphalow,k))/rwdsavings; % MT proportion of time savings
    rwdpropsyoung(2,k) = (alphacurvesyoung_rt(alphahigh,k) - alphacurvesyoung_rt(alphalow,k))/rwdsavings; % RT proportion of time savings
    
    alphacurvesold=squeeze(old(:,effscaleind,k,probscaleind,:));
    alphacurvesold_mt(:,k)=alphacurvesold(:,1);
    alphacurvesold_rt(:,k)=alphacurvesold(:,2);
    rwdsavingsold = ((alphacurvesold_rt(alphahigh,k) - alphacurvesold_rt(alphalow,k))+(alphacurvesold_mt(alphahigh,k) - alphacurvesold_mt(alphalow,k)));
    rwddeltasold(1,k) = (alphacurvesold_mt(alphahigh,k) - alphacurvesold_mt(alphalow,k));
    rwddeltasold(2,k) = (alphacurvesold_rt(alphahigh,k) - alphacurvesold_rt(alphalow,k));
    rwdpropsold(1,k) = (alphacurvesold_mt(alphahigh,k) - alphacurvesold_mt(alphalow,k))/rwdsavingsold; % MT proportion of time savings
    rwdpropsold(2,k) = (alphacurvesold_rt(alphahigh,k) - alphacurvesold_rt(alphalow,k))/rwdsavingsold; % RT proportion of time savings
end

% probability scale differences
for m=1:length(myprobscales)
    probcurvesyoung=squeeze(young(:,effscaleind,alphascaleind,m,:));
    probcurvesyoung_mt(:,m)=probcurvesyoung(:,1);
    probcurvesyoung_rt(:,m)=probcurvesyoung(:,2);
    probsavings = ((probcurvesyoung_rt(alphahigh,m) - probcurvesyoung_rt(alphalow,m))+(probcurvesyoung_mt(alphahigh,m) - probcurvesyoung_mt(alphalow,m)));
    probdeltasyoung(1,m) = (probcurvesyoung_mt(alphahigh,m) - probcurvesyoung_mt(alphalow,m));
    probdeltasyoung(2,m) = (probcurvesyoung_rt(alphahigh,m) - probcurvesyoung_rt(alphalow,m));
    probpropsyoung(1,m) = (probcurvesyoung_mt(alphahigh,m) - probcurvesyoung_mt(alphalow,m))/probsavings; % MT proportion of time savings
    probpropsyoung(2,m) = (probcurvesyoung_rt(alphahigh,m) - probcurvesyoung_rt(alphalow,m))/probsavings; % RT proportion of time savings
    
    probcurvesold=squeeze(old(:,effscaleind,alphascaleind,m,:));
    probcurvesold_mt(:,m)=probcurvesold(:,1);
    probcurvesold_rt(:,m)=probcurvesold(:,2);
    probsavingsold = ((probcurvesold_rt(alphahigh,m) - probcurvesold_rt(alphalow,m))+(probcurvesold_mt(alphahigh,m) - probcurvesold_mt(alphalow,m)));
    probdeltasold(1,m) = (probcurvesold_mt(alphahigh,m) - probcurvesold_mt(alphalow,m));
    probdeltasold(2,m) = (probcurvesold_rt(alphahigh,m) - probcurvesold_rt(alphalow,m));
    probpropsold(1,m) = (probcurvesold_mt(alphahigh,m) - probcurvesold_mt(alphalow,m))/probsavingsold; % MT proportion of time savings
    probpropsold(2,m) = (probcurvesold_rt(alphahigh,m) - probcurvesold_rt(alphalow,m))/probsavingsold; % RT proportion of time savings
end

% Proportions of savings from empirical data
dattimes(:,1) = table2array(aggdat(aggdat.MASS ==0 & aggdat.REWD == 1, [14,16]))' - table2array(aggdat(aggdat.MASS ==0 & aggdat.REWD == 0, [14,16]))'; % col 1 = no mass, row 1 = deltatRT
dattimes(:,2) = table2array(aggdat(aggdat.MASS ==1 & aggdat.REWD == 1, [14,16]))' - table2array(aggdat(aggdat.MASS ==1 & aggdat.REWD == 0, [14,16]))';
dattimes(:,3) = [0.286 - 0.3015; 0.567 - 0.5946]; % young adults from Erik
dattimes(:,4) = [0.396 - 0.417; 0.672 - 0.687]; % old adults from Erik

datprops(1,1) = dattimes(2,1)/sum(dattimes(:,1)); datprops(2,1) = 1-datprops(1,1);
datprops(1,2) = dattimes(2,2)/sum(dattimes(:,2)); datprops(2,2) = 1-datprops(1,2);
datprops(1,3) = dattimes(2,3)/sum(dattimes(:,3)); datprops(2,3) = 1-datprops(1,3);
datprops(1,4) = dattimes(2,4)/sum(dattimes(:,4)); datprops(2,4) = 1-datprops(1,4);


% For 3D plots
RT = zeros(5);
MT = zeros(5);
RTprops = zeros(5);
totalsavings = zeros(5);

for k=1:length(myalphascales)
    for j=1:length(myeffscales)
        % young data
        curvesyoung=squeeze(young(:,j,k,probscaleind,:));
        curvesyoung_mt(:,j)=curvesyoung(:,1);
        curvesyoung_rt(:,j)=curvesyoung(:,2);
        
        savings = (curvesyoung_rt(alphahigh,j) - curvesyoung_rt(alphalow,j)) + (curvesyoung_mt(alphahigh,j) - curvesyoung_mt(alphalow,j)); % total time savings
        RT(j,k) = (curvesyoung_rt(alphahigh,j) - curvesyoung_rt(alphalow,j));
        MT(j,k) = (curvesyoung_mt(alphahigh,j) - curvesyoung_mt(alphalow,j));
        RTprops(j,k) = (curvesyoung_rt(alphahigh,j) - curvesyoung_rt(alphalow,j))/savings;
        totalsavings(j,k) = savings; %./(curvesyoung_rt(alphalow,j)+curvesyoung_mt(alphalow,j));
    end
end

effscalemat = repmat(myeffscales',1,5);
alphascalemat = repmat(myalphascales,5,1);
C = (1./effscalemat).*alphascalemat;

%% Figure for time savings, stacked
figure
%     subplot(4,1,1)
%         [tme, r, v, a] = minjerk([0,0],[0,0.1],1,.001);
%         plot(tme,v(:,2))

    subplot(4,1,1)
        barh(categorical(myeffscales), [effdeltasyoung(2,:)' effdeltasyoung(1,:)'],'stacked');
        ylabel(sprintf('Increasing Effort'));
        xlabel(sprintf(['Time savings (s)\n%dJ - %dJ'], rwdhigh, rwdlow));
        xlim([-inf 0])
        set(gca,'ydir','reverse','xdir','reverse');
        legend('RT','MT','Location','southeast');
    subplot(4,1,2)
        barh(categorical(myalphascales), [rwddeltasyoung(2,:)' rwddeltasyoung(1,:)'] ,'stacked');
        ylabel(sprintf('Reducing Reward'));
        xlabel(sprintf(['Time savings (s)\n%dJ - %dJ'], rwdhigh, rwdlow));
        set(gca,'xdir','reverse')
    subplot(4,1,3)
        barh(categorical(myprobscales), [probdeltasyoung(2,:)' probdeltasyoung(1,:)'] ,'stacked');
        ylabel(sprintf('Reducing Accuracy'));
        xlabel(sprintf(['Time savings (s)\n%dJ - %dJ'], rwdhigh, rwdlow));
        set(gca,'ydir','reverse','xdir','reverse');
    subplot(4,1,4)
        x = categorical({'LOW','HIGH' 'YOUNG' 'OLD'});
        barh(reordercats(x, string(x)), [dattimes(2,:)' dattimes(1,:)'] ,'stacked');
        ylabel('Effort');
        xlabel(sprintf('Time savings (s)\nRWD - NRWD (a.u.)'));
        set(gca,'xdir','reverse')
    beautifyfig;
    
%% Figure for proportion of savings, stacked
figure
    subplot(4,1,1)
        barh(categorical(myeffscales), [effpropsyoung(2,:)' effpropsyoung(1,:)'],'stacked');
        ylabel(sprintf('Increasing Effort'));
        xlabel(sprintf(['Proportion of savings\n%dJ - %dJ'], rwdhigh, rwdlow));
        xlim([0 1])
        set(gca,'ydir','reverse');
        legend('RT','MT','Location','southeast');
    subplot(4,1,2)
        barh(categorical(myalphascales), [rwdpropsyoung(2,:)' rwdpropsyoung(1,:)'] ,'stacked');
        ylabel(sprintf('Reducing Reward'));
        xlabel(sprintf(['Proportion of savings\n%dJ - %dJ'], rwdhigh, rwdlow));
    subplot(4,1,3)
        barh(categorical(myprobscales), [probpropsyoung(2,:)' probpropsyoung(1,:)'] ,'stacked');
        ylabel(sprintf('Reducing Accuracy'));
        xlabel(sprintf(['Proportion of savings\n%dJ - %dJ'], rwdhigh, rwdlow));
        set(gca,'ydir','reverse');
    subplot(4,1,4)
        x = categorical({'LOW','HIGH' 'YOUNG' 'OLD'});
        barh(reordercats(x, string(x)), [datprops(2,:)' datprops(1,:)'] ,'stacked');
        ylabel('Effort');
        xlabel(sprintf('Proportion of savings\nRWD - NRWD (a.u.)'));
    beautifyfig;

%% RT/MT Speed Accuracy, Reward Curves, and Proportions

% Speed-Accuracy Stuff for Plotting
    % Scaled up RT probability function
    Pr_oldyoung = @(Tr) (1./(1 + exp(-(param.myc0to) - (param.myc1to).*Tr)));
    % Alternative RT function
    k00 = 20; 
    x00 = 0.30;
    Pr_test2 = @(Tr) (0.75./(1 + exp(-k00.*(Tr - x00))))+0.25;

    % Movement time fits
    load speedaccdata.mat
    % young fit
    [logitCoefyoung,devyoung] = glmfit(mt.young(:,1),pR.young.y(:,1),'binomial','logit');
    logitFityoung = glmval(logitCoefyoung,mt.young(:,1),'logit');
    % old fit
    [logitCoefold,devold] = glmfit(mt.old(:,1),pR.old.y(:,1),'binomial','logit');
    logitFitold = glmval(logitCoefold,mt.old(:,1),'logit');

% for shaded/highlighted patch of plot
    ptchidx = (myalphas >= rwdlow) & (myalphas <= rwdhigh);
    ptchcol = [221,160,221]./256;
    colgray = [0.6    0.6    0.6];
    unit = 1000; %scaling factor (to go from s to ms)

figure
clear h    
tiledlayout(7,2)

% RT speed-accuracy logisitc/CDFs
nexttile(1,[3 1])
    fplot(Pr_oldyoung,'Color','k', 'LineWidth', 1, 'DisplayName', 'Orig. RT function')
    hold on
    fplot(Pr_test2, 'Color', [.7 .7 .7], 'LineWidth', 1, 'DisplayName', 'Alt. RT function')
    xline(0.150,'k--','HandleVisibility','off'); yline(.25,'k--','HandleVisibility','off');
    legend('show', 'Location','southeast');
    xlim([-0.5 1.0]);
    xlabel('Reaction time (s)'); ylabel('Probability of success');
    hold off

% Fitted MT speed accuracy logisitc/CDFs    
nexttile(2,[3 1])
    scatter(mt.young(:,1),pR.young.y(:,1),'gs', 'MarkerFaceColor','g', 'MarkerFaceAlpha',0.5);
    hold on
    scatter(mt.old(:,1),pR.old.y(:,1),'rs','MarkerFaceColor','r', 'MarkerFaceAlpha',0.5);
    %plot models fits
    mts=0:0.05:1.5;
    logitCoefyoung =[ -4.9277 9.8371];
    logitCoefold =[ -8.8805 14.5769];
    % young fit
    b0=logitCoefyoung(1); b1=logitCoefyoung(2);
    plot(mts, 1./(1+exp(-b0 - mts*b1)),'g','LineWidth', 1);
    % old fit
    b0=logitCoefold(1); b1=logitCoefold(2);
    plot(mts, 1./(1+exp(-b0 - mts*b1)),'r','LineWidth', 1);
    yline(.5,'k--','HandleVisibility','off');
    hold off
    xlabel('Movement time (s)'); 
    ylabel('Probability of success');
    legend('Young','Old','Location','southeast');
    set(gca,'xlim',[0 1.5])

% RT/MT curves for effort scale
nexttile(7,[3 1])
    % young effort curves   
    h(:,1)=plot(myalphas,squeeze(young(:,3,alphascaleind,probscaleind,:)).*unit,'LineWidth', 1);
    hold on
    h(2,1).LineStyle = '--'; % dashed = RT
    h(1,1).Color = 'g';
    h(2,1).Color = 'g';
    % old effort curves
    h(:,2)=plot(myalphas,squeeze(old(:,3,alphascaleind,probscaleind,:)).*unit,'Color', 'r', 'LineWidth',1);
    h(2,2).LineStyle = '--';
    h(1,2).Color = 'r';
    h(2,2).Color = 'r';
    title('Increasing reaching effort')
    plot(repmat(myalphas(alphalow),[4 1]), [squeeze(old(alphalow,3,alphascaleind,probscaleind,:)); squeeze(young(alphalow,3,alphascaleind,probscaleind,:))].*unit, 'k>','MarkerSize',4.0)
    plot(repmat(myalphas(alphahigh),[4 1]), [squeeze(old(alphahigh,3,alphascaleind,probscaleind,:)); squeeze(young(alphahigh,3,alphascaleind,probscaleind,:))].*unit, 'k<','MarkerSize',4.0)
    % highlighted region
    patch([myalphas(ptchidx) fliplr(myalphas(ptchidx))], [1.5.*ones(size(myalphas(ptchidx))).*unit -0.5.*ones(size(myalphas(ptchidx))).*unit],...
    ptchcol, 'FaceAlpha',0.3, 'EdgeColor','none')    
    legend(h(:,1),'MT','RT');
    xlabel('Reward (alpha,J)')
    xticks([40 60 80 100])
    yticks([300 600 900 1200 1500])
    ylabel('Duration (ms)')
    set(gca,'ylim',[0.3 1.5].*unit,'xlim',[25 65])

% RT/MT curves for reward scale    
nexttile(8,[3 1])
    % 1.0 RWD scale curves   
    h(:,1)=plot(myalphas,squeeze(young(:,effscaleind,3,probscaleind,:)).*unit,'LineWidth', 1);
    hold on
    h(2,1).LineStyle = '--'; % dashed = RT
    h(1,1).Color = 'g';
    h(2,1).Color = 'g';
    % 0.8 RWD scale effort curves
    h(:,2)=plot(myalphas,squeeze(young(:,effscaleind,1,probscaleind,:)).*unit,...
        'Color', 'r', 'LineWidth',1);
    h(2,2).LineStyle = '--';
    h(1,2).Color = 'r';
    h(2,2).Color = 'r';
    title('Reducing reward scaling')
    plot(repmat(myalphas(alphalow),[4 1]), [squeeze(young(alphalow,effscaleind,3,probscaleind,:)); squeeze(young(alphalow,effscaleind,1,probscaleind,:))].*unit, 'k>','MarkerSize',4.0)
    plot(repmat(myalphas(alphahigh),[4 1]), [squeeze(young(alphahigh,effscaleind,3,probscaleind,:)); squeeze(young(alphahigh,effscaleind,1,probscaleind,:))].*unit, 'k<','MarkerSize',4.0)
    % highlighted region
    patch([myalphas(ptchidx) fliplr(myalphas(ptchidx))], [1.5.*ones(size(myalphas(ptchidx))).*unit -0.5.*ones(size(myalphas(ptchidx))).*unit],...
    ptchcol, 'FaceAlpha', 0.3, 'EdgeColor','none')
    legend(h(:,1),'MT','RT');
    xlabel('Reward (alpha,J)')
    xticks([40 60 80 100])
    yticks([300 600 900 1200 1500])
    ylabel('Duration (ms)')
    set(gca,'ylim',[0.3 1.5].*unit,'xlim',[25 65])

% Proportion savings barplots for effort scale
nexttile(13)
    b = barh(categorical({'Young' 'Old'}), [effpropsyoung(2,3)' effpropsyoung(1,3)'; effpropsold(2,3)' effpropsold(1,3)'],'stacked');
    b(1).FaceColor = 'flat';
    b(1).CData = [2/5 0 0; 0 2/5 0;];
    b(2).FaceColor = 'flat';
    b(2).CData = [5/5 0 0; 0 5/5 0];
    ylabel(sprintf('Increased effort'));
    xlabel(sprintf(['Proportion of savings\n%dJ - %dJ'], rwdhigh, rwdlow));
    
% Proportion savings barplots for reward scale
nexttile(14)
    b = barh(categorical({'High' 'Low'}), [rwdpropsyoung(2,3)' rwdpropsyoung(1,3)'; rwdpropsyoung(2,1)' rwdpropsyoung(1,1)'],'stacked');
    b(1).FaceColor = 'flat';
    b(1).CData = [0 2/5 0; 2/5 0 0];
    b(2).FaceColor = 'flat';
    b(2).CData = [0 5/5 0; 5/5 0 0];
    ylabel(sprintf('Reward scale'));
    xlabel(sprintf(['Proportion of savings\n%dJ - %dJ'], rwdhigh, rwdlow));
    set(gca,'ydir','reverse')

beautifyfig;

%% Combined figure using old/young data and 3D Surface plots

figure
clear h
% ADD 50-55J SHADING BOX
% 2x2 for first figure, 2x1 for RT/MT bars, 2x2 for last line graphs

% for shaded/highlighted patch of plot
ptchidx = (myalphas >= rwdlow) & (myalphas <= rwdhigh);
% ptchcol = [221,160,221]./256;
ptchcol = spring(1);
colgray = [0.6    0.6    0.6];
unit = 1000; %scaling factor (to go from s to ms)

tiledlayout(6,6)
nexttile(1,[4 3])
    for i = 1:5
        h(:,i)=plot(myalphas,squeeze(young(:,i,alphascaleind,probscaleind,:)).*unit,'LineWidth', 1);
        hold on
        h(1,i).Color = 1/i.*[0 1 0];
        h(2,i).Color = 1/i.*[0 1 0];
        h(2,i).LineStyle = '--';
    end
%     h(:,i+1)=plot(myalphas,squeeze(old(:,3,alphascaleind,probscaleind,:)).*unit,...
%         'Color', 'r', 'LineWidth',1);
%     h(2,i+1).LineStyle = '--';
%     legend(h(1,[1 3 5 6]),'Young 0.8', 'Young 1.0', 'Young 1.2', 'Old')
    legend(h(1,[1 3 5]),'Effort 0.8', 'Effort 1.0', 'Effort 1.2')
    
    patch([myalphas(ptchidx) fliplr(myalphas(ptchidx))], [1.5.*ones(size(myalphas(ptchidx))).*unit -0.5.*ones(size(myalphas(ptchidx))).*unit],...
    ptchcol, 'FaceAlpha',0.3, 'EdgeColor','none','HandleVisibility','off')
    xlabel('Reward (alpha,J)')
    xticks([40 60 80 100])
    yticks([300 600 900 1200 1500])
    ylabel('Duration (ms)')
    set(gca,'ylim',[0.3 1.2].*unit,'xlim',[28 100])
    hold off
    
nexttile(4,[4,3])
    for i = 1:5
        h(:,i)=plot(myalphas,squeeze(young(:,effscaleind,i,probscaleind,:)).*unit,...
                    'LineWidth', 1);
        hold on
        h(2,i).LineStyle = '--';
        h(1,i).Color = i.*[0 1/5 0];
        h(2,i).Color = i.*[0 1/5 0];
    end
    
    patch([myalphas(ptchidx) fliplr(myalphas(ptchidx))], [1.5.*ones(size(myalphas(ptchidx))).*unit -0.5.*ones(size(myalphas(ptchidx))).*unit],...
        ptchcol, 'FaceAlpha',0.3, 'EdgeColor','none')
    legend(h(1,[1 3 5]),'RWD 0.8', 'RWD 1.0', 'RWD 1.2')
    xlabel('Reward (alpha,J)')
    ylabel('Duration (ms)')
    xticks([40 60 80 100])
    yticks([300 600 900 1200 1500])
    set(gca,'ylim',[0.3 1.2].*unit,'xlim',[28 100])
    hold off
    
nexttile(25,[2,2])
    surf(effscalemat, alphascalemat, RT, C, 'FaceAlpha', 0.5);
    xlabel('effort scale'); ylabel('reward scale'); zlabel('\DeltaRT');
    colormap('spring');

nexttile(27, [2,2])
    surf(effscalemat, alphascalemat, MT, C, 'FaceAlpha', 0.5);
    xlabel('effort scale'); ylabel('reward scale'); zlabel('\DeltaMT');
    colormap('spring');

nexttile(29,[2,2])
    surf(effscalemat, alphascalemat, RTprops, C, 'FaceAlpha', 0.5);
    zlim([0 1]);
    xlabel('effort scale'); ylabel('reward scale'); zlabel('RT proportion');
    colormap('spring');

cb = colorbar;
cb.Layout.Tile = 'south';
cb.TickDirection = 'out';
cb.Ticks = [min(C(:)), (min(C(:))+max(C(:)))/2, max(C(:))];
cb.TickLabels = {'Low' 'Med' 'High'};
cb.Label.String = 'Relative Utility';

beautifyfig;

%% RT/MT Speed Accuracy Changes

% for shaded/highlighted patch of plot
    ptchidx = (myalphas >= rwdlow) & (myalphas <= rwdhigh);
    ptchcol = [221,160,221]./256;
    colgray = [0.6    0.6    0.6];
    unit = 1000; %scaling factor (to go from s to ms)

figure
clear h    
tiledlayout(7,1)    
nexttile(1,[3 1])
    mts = 0:0.05:1.5;
    plot(mts, 1./(1+exp(-(-5) - mts*(10))),'g', 'LineWidth',1);
    hold on
    plot(mts, 1./(1+exp(-(-5) - mts*(7.5))),'r', 'LineWidth',1);
    % plot(mts, 1./(1+exp(-(-9) - mts*(15))),'g');
    % plot(mts, 1./(1+exp(-(-6) - mts*(8))),'g--');
    xlabel('Movement time (s)'); ylabel('P(success)'); title('Changing speed-accuracy tradeoff for MT')
    legend('High Acc.: c_0 = -5, c_1 = 10','Low Acc.: c_0 = -5, c_1 = 7.5', 'Location','southeast')

nexttile(4,[3 1])
    % young effort curves   
    h(:,1)=plot(myalphas,squeeze(young(:,3,alphascaleind,probscaleind,:)).*unit,'LineWidth', 1);
    hold on
    h(2,1).LineStyle = '--'; % dashed = RT
    h(1,1).Color = 'g';
    h(2,1).Color = 'g';
    % old effort curves
    h(:,2)=plot(myalphas,squeeze(old(:,3,alphascaleind,probscaleind,:)).*unit,'Color', 'r', 'LineWidth',1);
    h(2,2).LineStyle = '--';
    h(1,2).Color = 'r';
    h(2,2).Color = 'r';
    title('Increasing reaching effort')
    plot(repmat(myalphas(alphalow),[4 1]), [squeeze(old(alphalow,3,alphascaleind,probscaleind,:)); squeeze(young(alphalow,3,alphascaleind,probscaleind,:))].*unit, 'k>','MarkerSize',4.0)
    plot(repmat(myalphas(alphahigh),[4 1]), [squeeze(old(alphahigh,3,alphascaleind,probscaleind,:)); squeeze(young(alphahigh,3,alphascaleind,probscaleind,:))].*unit, 'k<','MarkerSize',4.0)
    % highlighted region
    patch([myalphas(ptchidx) fliplr(myalphas(ptchidx))], [1.5.*ones(size(myalphas(ptchidx))).*unit -0.5.*ones(size(myalphas(ptchidx))).*unit],...
    ptchcol, 'FaceAlpha',0.3, 'EdgeColor','none')    
    legend(h(:,1),'MT','RT');
    xlabel('Reward (alpha,J)')
    xticks([40 60 80 100])
%     yticks([300 600 900 1200 1500])
    ylabel('Duration (ms)')
%     set(gca,'ylim',[0.3 1.5].*unit,'xlim',[25 65])
    set(gca,'xlim',[25 65])

        
nexttile(7)
    clear b
    b = barh(categorical({'High Acc.','Low Acc.'}), [effpropsyoung(2,3) effpropsyoung(1,3); effpropsold(2,3) effpropsold(1,3)], 'stacked');
    b(1).FaceColor = 'flat';
    b(1).CData = [0 2/5 0; 2/5 0 0];
    b(2).FaceColor = 'flat';
    b(2).CData = [0 1 0; 1 0 0];
    ylabel(sprintf('Increased effort'));
    xlabel(sprintf(['Proportion of savings\n%dJ - %dJ'], rwdhigh, rwdlow));

beautifyfig;

%% Effort, Reward, and Probability Scaling combined

% Speed-Accuracy Stuff for Plotting
    % Scaled up RT probability function
    Pr_oldyoung = @(Tr) (1./(1 + exp(-(param.myc0to) - (param.myc1to).*Tr)));
    % Alternative RT function
    k00 = 20; 
    x00 = 0.30;
    Pr_test2 = @(Tr) (0.75./(1 + exp(-k00.*(Tr - x00))))+0.25;

    % Movement time fits
    load speedaccdata.mat
    % young fit
    [logitCoefyoung,devyoung] = glmfit(mt.young(:,1),pR.young.y(:,1),'binomial','logit');
    logitFityoung = glmval(logitCoefyoung,mt.young(:,1),'logit');
    % old fit
    [logitCoefold,devold] = glmfit(mt.old(:,1),pR.old.y(:,1),'binomial','logit');
    logitFitold = glmval(logitCoefold,mt.old(:,1),'logit');
    
    % Movement time speed-accuracy with changing probability scale
    logitAcc = @(Tr) 1./(1+exp(-(param.myc0) - Tr.*(param.myc1)));
    logitInacc = @(Tr) 1./(1+exp(-(1.2.*param.myc0) - Tr.*(param.myc1))); % higher probability scale = more inaccuracy
    
% for shaded/highlighted patch of plot
    ptchidx = (myalphas >= rwdlow) & (myalphas <= rwdhigh);
    ptchcol = [221,160,221]./256;
    colgray = [0.6    0.6    0.6];
    unit = 1000; %scaling factor (to go from s to ms)

figure
clear h    
tiledlayout(8,3)

% RT speed-accuracy logisitic/CDFs
nexttile(1,[3 1])
    fplot(Pr_oldyoung,'Color','k', 'LineWidth', 1, 'DisplayName', 'Orig. RT function')
    hold on
    fplot(Pr_test2, 'Color', [.7 .7 .7], 'LineWidth', 1, 'DisplayName', 'Alt. RT function')
    xline(0.150,'k--','HandleVisibility','off'); yline(.25,'k--','HandleVisibility','off');
    legend('show', 'Location','southeast');
    xlim([-0.5 1.0]);
    xlabel('Reaction time (s)'); ylabel('Probability of success');
    hold off

% Fitted MT speed-accuracy logisitic curves
nexttile(2,[3 1])
    scatter(mt.young(:,1),pR.young.y(:,1),'gs', 'MarkerFaceColor','g', 'MarkerFaceAlpha',0.5);
    hold on
    scatter(mt.old(:,1),pR.old.y(:,1),'rs','MarkerFaceColor','r', 'MarkerFaceAlpha',0.5);
    %plot models fits
    mts=0:0.05:1.5;
    logitCoefyoung =[ -4.9277 9.8371];
    logitCoefold =[ -8.8805 14.5769];
    % young fit
    b0=logitCoefyoung(1); b1=logitCoefyoung(2);
    plot(mts, 1./(1+exp(-b0 - mts*b1)),'g','LineWidth', 1);
    % old fit
    b0=logitCoefold(1); b1=logitCoefold(2);
    plot(mts, 1./(1+exp(-b0 - mts*b1)),'r','LineWidth', 1);
    yline(.5,'k--','HandleVisibility','off');
    hold off
    xlabel('Movement time (s)'); 
    ylabel('Probability of success');
    legend('Young','Old','Location','southeast');
    set(gca,'xlim',[0 1.5])

% MT speed-accuracy curves for changing probability scale
nexttile(3,[3 1])
    fplot(logitAcc,'Color','g', 'LineWidth', 1, 'DisplayName', 'Prob. Scale = 1.0')
    hold on
    fplot(logitInacc, 'Color', 'r', 'LineWidth', 1, 'DisplayName', 'Prob. Scale = 1.2')
    legend('show', 'Location','southeast');
    xlim([0 1.5]);
    xlabel('Movement time (s)'); ylabel('Probability of success');
    hold off

% RT/MT curves for effort scale
nexttile(10,[3 1])
    % young effort curves   
    h(:,1)=plot(myalphas,squeeze(young(:,3,alphascaleind,probscaleind,:)).*unit,'LineWidth', 1);
    hold on
    h(2,1).LineStyle = '--'; % dashed = RT
    h(1,1).Color = 'g';
    h(2,1).Color = 'g';
    % old effort curves
    h(:,2)=plot(myalphas,squeeze(old(:,3,alphascaleind,probscaleind,:)).*unit,'Color', 'r', 'LineWidth',1);
    h(2,2).LineStyle = '--';
    h(1,2).Color = 'r';
    h(2,2).Color = 'r';
    title('Increasing reaching effort')
    plot(repmat(myalphas(alphalow),[4 1]), [squeeze(old(alphalow,3,alphascaleind,probscaleind,:)); squeeze(young(alphalow,3,alphascaleind,probscaleind,:))].*unit, 'k>','MarkerSize',4.0)
    plot(repmat(myalphas(alphahigh),[4 1]), [squeeze(old(alphahigh,3,alphascaleind,probscaleind,:)); squeeze(young(alphahigh,3,alphascaleind,probscaleind,:))].*unit, 'k<','MarkerSize',4.0)
    % highlighted region
    patch([myalphas(ptchidx) fliplr(myalphas(ptchidx))], [1.5.*ones(size(myalphas(ptchidx))).*unit -0.5.*ones(size(myalphas(ptchidx))).*unit],...
    ptchcol, 'FaceAlpha',0.3, 'EdgeColor','none')    
    legend(h(:,1),'MT','RT');
    xlabel('Reward (alpha,J)')
    xticks([40 60 80 100])
    yticks([300 600 900 1200 1500])
    ylabel('Duration (ms)')
    set(gca,'ylim',[0.3 1.5].*unit,'xlim',[25 65])
    
% RT/MT curves for reward scale
nexttile(11,[3 1])
    % 1.0 RWD scale curves   
    h(:,1)=plot(myalphas,squeeze(young(:,effscaleind,3,probscaleind,:)).*unit,'LineWidth', 1);
    hold on
    h(2,1).LineStyle = '--'; % dashed = RT
    h(1,1).Color = 'g';
    h(2,1).Color = 'g';
    % 0.8 RWD scale effort curves
    h(:,2)=plot(myalphas,squeeze(young(:,effscaleind,1,probscaleind,:)).*unit,...
        'Color', 'r', 'LineWidth',1);
    h(2,2).LineStyle = '--';
    h(1,2).Color = 'r';
    h(2,2).Color = 'r';
    title('Reducing reward scaling')
    plot(repmat(myalphas(alphalow),[4 1]), [squeeze(young(alphalow,effscaleind,3,probscaleind,:)); squeeze(young(alphalow,effscaleind,1,probscaleind,:))].*unit, 'k>','MarkerSize',4.0)
    plot(repmat(myalphas(alphahigh),[4 1]), [squeeze(young(alphahigh,effscaleind,3,probscaleind,:)); squeeze(young(alphahigh,effscaleind,1,probscaleind,:))].*unit, 'k<','MarkerSize',4.0)
    % highlighted region
    patch([myalphas(ptchidx) fliplr(myalphas(ptchidx))], [1.5.*ones(size(myalphas(ptchidx))).*unit -0.5.*ones(size(myalphas(ptchidx))).*unit],...
    ptchcol, 'FaceAlpha', 0.3, 'EdgeColor','none')
    legend(h(:,1),'MT','RT');
    xlabel('Reward (alpha,J)')
    xticks([40 60 80 100])
    yticks([300 600 900 1200 1500])
    ylabel('Duration (ms)')
    set(gca,'ylim',[0.3 1.5].*unit,'xlim',[25 65])

% RT/MT curves for probability scale
nexttile(12,[3 1])
    % 1.0 Prob scale curves   
    h(:,1)=plot(myalphas,squeeze(young(:,effscaleind,alphascaleind,3,:)).*unit,'LineWidth', 1);
    hold on
    h(2,1).LineStyle = '--'; % dashed = RT
    h(1,1).Color = 'g';
    h(2,1).Color = 'g';
    % 1.2 Prob scale effort curves
    h(:,2)=plot(myalphas,squeeze(young(:,effscaleind,alphascaleind,5,:)).*unit,...
        'Color', 'r', 'LineWidth',1);
    h(2,2).LineStyle = '--';
    h(1,2).Color = 'r';
    h(2,2).Color = 'r';
    title('Reducing accuracy')
    plot(repmat(myalphas(alphalow),[4 1]), [squeeze(young(alphalow,effscaleind,alphascaleind,3,:)); squeeze(young(alphalow,effscaleind,alphascaleind,5,:))].*unit, 'k>','MarkerSize',4.0)
    plot(repmat(myalphas(alphahigh),[4 1]), [squeeze(young(alphahigh,effscaleind,alphascaleind,3,:)); squeeze(young(alphahigh,effscaleind,alphascaleind,5,:))].*unit, 'k<','MarkerSize',4.0)
    % highlighted region
    patch([myalphas(ptchidx) fliplr(myalphas(ptchidx))], [1.5.*ones(size(myalphas(ptchidx))).*unit -0.5.*ones(size(myalphas(ptchidx))).*unit],...
    ptchcol, 'FaceAlpha', 0.3, 'EdgeColor','none')
    legend(h(:,1),'MT','RT');
    xlabel('Reward (alpha,J)')
    xticks([40 60 80 100])
    yticks([300 600 900 1200 1500])
    ylabel('Duration (ms)')
    set(gca,'ylim',[0.3 1.5].*unit,'xlim',[25 65])
    
% nexttile(19)
%     b = barh(categorical({'Young' 'Old'}), [effpropsyoung(2,3)' effpropsyoung(1,3)'; effpropsold(2,3)' effpropsold(1,3)'],'stacked');
%     b(1).FaceColor = 'flat';
%     b(1).CData = [2/5 0 0; 0 2/5 0;];
%     b(2).FaceColor = 'flat';
%     b(2).CData = [5/5 0 0; 0 5/5 0];
%     ylabel(sprintf('Increased effort'));
%     xlabel(sprintf(['Proportion of savings\n%dJ - %dJ'], rwdhigh, rwdlow));
%     
% 
% nexttile(20)
%     b = barh(categorical({'High' 'Low'}), [rwdpropsyoung(2,3)' rwdpropsyoung(1,3)'; rwdpropsyoung(2,1)' rwdpropsyoung(1,1)'],'stacked');
%     b(1).FaceColor = 'flat';
%     b(1).CData = [0 2/5 0; 2/5 0 0];
%     b(2).FaceColor = 'flat';
%     b(2).CData = [0 5/5 0; 5/5 0 0];
%     ylabel(sprintf('Reward scale'));
%     xlabel(sprintf(['Proportion of savings\n%dJ - %dJ'], rwdhigh, rwdlow));
%     set(gca,'ydir','reverse')
%     
% nexttile(21)
%     b = barh(categorical({'High' 'Low'}), [probpropsyoung(2,3)' probpropsyoung(1,3)'; probpropsyoung(2,5)' probpropsyoung(1,5)'],'stacked');
%     b(1).FaceColor = 'flat';
%     b(1).CData = [0 2/5 0; 2/5 0 0];
%     b(2).FaceColor = 'flat';
%     b(2).CData = [0 5/5 0; 5/5 0 0];
%     ylabel(sprintf('Accuracy scale'));
%     xlabel(sprintf(['Proportion of savings\n%dJ - %dJ'], rwdhigh, rwdlow));
%     set(gca,'ydir','reverse')

% Time savings barplots for effort scale
nexttile(19)
    b = barh(categorical({'Young' 'Old'}), [effdeltasyoung(2,3)' effdeltasyoung(1,3)'; effdeltasold(2,3)' effdeltasold(1,3)'],'stacked');
    b(1).FaceColor = 'flat';
    b(1).CData = [2/5 0 0; 0 2/5 0;];
    b(2).FaceColor = 'flat';
    b(2).CData = [5/5 0 0; 0 5/5 0];
    ylabel(sprintf('Increased effort'));
    xlabel(sprintf(['Time savings (s)\n%dJ - %dJ'], rwdhigh, rwdlow));
    set(gca, 'xdir','reverse')

% Time savings barplots for reward scale 
nexttile(20)
    b = barh(categorical({'High' 'Low'}), [rwddeltasyoung(2,3)' rwddeltasyoung(1,3)'; rwddeltasyoung(2,1)' rwddeltasyoung(1,1)'],'stacked');
    b(1).FaceColor = 'flat';
    b(1).CData = [0 2/5 0; 2/5 0 0];
    b(2).FaceColor = 'flat';
    b(2).CData = [0 5/5 0; 5/5 0 0];
    ylabel(sprintf('Reward scale'));
    xlabel(sprintf(['Time savings (s)\n%dJ - %dJ'], rwdhigh, rwdlow));
    set(gca,'ydir','reverse', 'xdir','reverse')   

% Time savings barplots for probability scale
nexttile(21)
    b = barh(categorical({'High' 'Low'}), [probdeltasyoung(2,3)' probdeltasyoung(1,3)'; probdeltasyoung(2,5)' probdeltasyoung(1,5)'],'stacked');
    b(1).FaceColor = 'flat';
    b(1).CData = [0 2/5 0; 2/5 0 0];
    b(2).FaceColor = 'flat';
    b(2).CData = [0 5/5 0; 5/5 0 0];
    ylabel(sprintf('Accuracy scale'));
    xlabel(sprintf(['Time savings (s)\n%dJ - %dJ'], rwdhigh, rwdlow));
    set(gca,'ydir','reverse', 'xdir','reverse') 

% Proportion time savings barplots for effort scale
nexttile(22)
    b = barh(categorical({'Young' 'Old'}), [effpropsyoung(2,3)' effpropsyoung(1,3)'; effpropsold(2,3)' effpropsold(1,3)'],'stacked');
    b(1).FaceColor = 'flat';
    b(1).CData = [2/5 0 0; 0 2/5 0;];
    b(2).FaceColor = 'flat';
    b(2).CData = [5/5 0 0; 0 5/5 0];
    ylabel(sprintf('Increased effort'));
    xlabel(sprintf(['Proportion of time savings\n%dJ - %dJ'], rwdhigh, rwdlow));

% Proportion time savings barplots for reward scale 
nexttile(23)
    b = barh(categorical({'High' 'Low'}), [rwdpropsyoung(2,3)' rwdpropsyoung(1,3)'; rwdpropsyoung(2,1)' rwdpropsyoung(1,1)'],'stacked');
    b(1).FaceColor = 'flat';
    b(1).CData = [0 2/5 0; 2/5 0 0];
    b(2).FaceColor = 'flat';
    b(2).CData = [0 5/5 0; 5/5 0 0];
    ylabel(sprintf('Reward scale'));
    xlabel(sprintf(['Proportion of time savings\n%dJ - %dJ'], rwdhigh, rwdlow));
    %set(gca,'ydir','reverse', 'xdir','reverse')   
    set(gca,'ydir','reverse') 

% Proportion time savings barplots for probability scale
nexttile(24)
    b = barh(categorical({'High' 'Low'}), [probpropsyoung(2,3)' probpropsyoung(1,3)'; probpropsyoung(2,5)' probpropsyoung(1,5)'],'stacked');
    b(1).FaceColor = 'flat';
    b(1).CData = [0 2/5 0; 2/5 0 0];
    b(2).FaceColor = 'flat';
    b(2).CData = [0 5/5 0; 5/5 0 0];
    ylabel(sprintf('Accuracy scale'));
    xlabel(sprintf(['Proportion of time savings\n%dJ - %dJ'], rwdhigh, rwdlow));
    %set(gca,'ydir','reverse', 'xdir','reverse') 
    set(gca,'ydir','reverse')

beautifyfig;

%% Effort, Reward, and Probability Scaling combined v2
% includes both delta savings and proportion savings

figure
clear h    
tiledlayout(9,3)    

% RT/MT curves for effort scale
nexttile(1,[3 1])
    % young effort curves   
    h(:,1)=plot(myalphas,squeeze(young(:,3,alphascaleind,probscaleind,:)).*unit,'LineWidth', 1);
    hold on
    h(2,1).LineStyle = '--'; % dashed = RT
    h(1,1).Color = 'g';
    h(2,1).Color = 'g';
    % old effort curves
    h(:,2)=plot(myalphas,squeeze(old(:,3,alphascaleind,probscaleind,:)).*unit,'Color', 'r', 'LineWidth',1);
    h(2,2).LineStyle = '--';
    h(1,2).Color = 'r';
    h(2,2).Color = 'r';
    title('Increasing reaching effort')
    plot(repmat(myalphas(alphalow),[4 1]), [squeeze(old(alphalow,3,alphascaleind,probscaleind,:)); squeeze(young(alphalow,3,alphascaleind,probscaleind,:))].*unit, 'k>','MarkerSize',4.0)
    plot(repmat(myalphas(alphahigh),[4 1]), [squeeze(old(alphahigh,3,alphascaleind,probscaleind,:)); squeeze(young(alphahigh,3,alphascaleind,probscaleind,:))].*unit, 'k<','MarkerSize',4.0)
    % highlighted region
    patch([myalphas(ptchidx) fliplr(myalphas(ptchidx))], [1.5.*ones(size(myalphas(ptchidx))).*unit -0.5.*ones(size(myalphas(ptchidx))).*unit],...
    ptchcol, 'FaceAlpha',0.3, 'EdgeColor','none')    
    legend(h(:,1),'MT','RT');
    xlabel('Reward (alpha,J)')
    xticks([40 60 80 100])
    yticks([300 600 900 1200 1500])
    ylabel('Duration (ms)')
    set(gca,'ylim',[0.3 1.5].*unit,'xlim',[25 65])

% RT/MT curves for reward scale
nexttile(10,[3 1])
    % 1.0 RWD scale curves   
    h(:,1)=plot(myalphas,squeeze(young(:,effscaleind,3,probscaleind,:)).*unit,'LineWidth', 1);
    hold on
    h(2,1).LineStyle = '--'; % dashed = RT
    h(1,1).Color = 'g';
    h(2,1).Color = 'g';
    % 0.8 RWD scale effort curves
    h(:,2)=plot(myalphas,squeeze(young(:,effscaleind,1,probscaleind,:)).*unit,...
        'Color', 'r', 'LineWidth',1);
    h(2,2).LineStyle = '--';
    h(1,2).Color = 'r';
    h(2,2).Color = 'r';
    title('Reducing reward scaling')
    plot(repmat(myalphas(alphalow),[4 1]), [squeeze(young(alphalow,effscaleind,3,probscaleind,:)); squeeze(young(alphalow,effscaleind,1,probscaleind,:))].*unit, 'k>','MarkerSize',4.0)
    plot(repmat(myalphas(alphahigh),[4 1]), [squeeze(young(alphahigh,effscaleind,3,probscaleind,:)); squeeze(young(alphahigh,effscaleind,1,probscaleind,:))].*unit, 'k<','MarkerSize',4.0)
    % highlighted region
    patch([myalphas(ptchidx) fliplr(myalphas(ptchidx))], [1.5.*ones(size(myalphas(ptchidx))).*unit -0.5.*ones(size(myalphas(ptchidx))).*unit],...
    ptchcol, 'FaceAlpha', 0.3, 'EdgeColor','none')
    legend(h(:,1),'MT','RT');
    xlabel('Reward (alpha,J)')
    xticks([40 60 80 100])
    yticks([300 600 900 1200 1500])
    ylabel('Duration (ms)')
    set(gca,'ylim',[0.3 1.5].*unit,'xlim',[25 65])

% RT/MT curves for probability scale    
nexttile(19,[3 1])
    % 1.0 Prob scale curves   
    h(:,1)=plot(myalphas,squeeze(young(:,effscaleind,alphascaleind,3,:)).*unit,'LineWidth', 1);
    hold on
    h(2,1).LineStyle = '--'; % dashed = RT
    h(1,1).Color = 'g';
    h(2,1).Color = 'g';
    % 1.2 Prob scale effort curves
    h(:,2)=plot(myalphas,squeeze(young(:,effscaleind,alphascaleind,5,:)).*unit,...
        'Color', 'r', 'LineWidth',1);
    h(2,2).LineStyle = '--';
    h(1,2).Color = 'r';
    h(2,2).Color = 'r';
    title('Reducing accuracy')
    plot(repmat(myalphas(alphalow),[4 1]), [squeeze(young(alphalow,effscaleind,alphascaleind,3,:)); squeeze(young(alphalow,effscaleind,alphascaleind,5,:))].*unit, 'k>','MarkerSize',4.0)
    plot(repmat(myalphas(alphahigh),[4 1]), [squeeze(young(alphahigh,effscaleind,alphascaleind,3,:)); squeeze(young(alphahigh,effscaleind,alphascaleind,5,:))].*unit, 'k<','MarkerSize',4.0)
    % highlighted region
    patch([myalphas(ptchidx) fliplr(myalphas(ptchidx))], [1.5.*ones(size(myalphas(ptchidx))).*unit -0.5.*ones(size(myalphas(ptchidx))).*unit],...
    ptchcol, 'FaceAlpha', 0.3, 'EdgeColor','none')
    legend(h(:,1),'MT','RT');
    xlabel('Reward (alpha,J)')
    xticks([40 60 80 100])
    yticks([300 600 900 1200 1500])
    ylabel('Duration (ms)')
    set(gca,'ylim',[0.3 1.5].*unit,'xlim',[25 65])

% Time savings barplots for effort scale    
nexttile(2, [3 1])
    b = barh(categorical({'Young' 'Old'}), [effdeltasyoung(2,3)' effdeltasyoung(1,3)'; effdeltasold(2,3)' effdeltasold(1,3)'],'stacked');
    b(1).FaceColor = 'flat';
    b(1).CData = [2/5 0 0; 0 2/5 0;];
    b(2).FaceColor = 'flat';
    b(2).CData = [5/5 0 0; 0 5/5 0];
    ylabel(sprintf('Increased effort'));
    xlabel(sprintf(['Time savings (s)\n%dJ - %dJ'], rwdhigh, rwdlow));
    set(gca, 'xdir','reverse')

% Time savings barplots for reward scale
nexttile(11, [3,1])
    b = barh(categorical({'High' 'Low'}), [rwddeltasyoung(2,3)' rwddeltasyoung(1,3)'; rwddeltasyoung(2,1)' rwddeltasyoung(1,1)'],'stacked');
    b(1).FaceColor = 'flat';
    b(1).CData = [0 2/5 0; 2/5 0 0];
    b(2).FaceColor = 'flat';
    b(2).CData = [0 5/5 0; 5/5 0 0];
    ylabel(sprintf('Reward scale'));
    xlabel(sprintf(['Time savings (s)\n%dJ - %dJ'], rwdhigh, rwdlow));
    set(gca,'ydir','reverse', 'xdir','reverse')

% Time savings barplots for probability scale
nexttile(20, [3,1])
    b = barh(categorical({'High' 'Low'}), [probdeltasyoung(2,3)' probdeltasyoung(1,3)'; probdeltasyoung(2,5)' probdeltasyoung(1,5)'],'stacked');
    b(1).FaceColor = 'flat';
    b(1).CData = [0 2/5 0; 2/5 0 0];
    b(2).FaceColor = 'flat';
    b(2).CData = [0 5/5 0; 5/5 0 0];
    ylabel(sprintf('Accuracy scale'));
    xlabel(sprintf(['Time savings (s)\n%dJ - %dJ'], rwdhigh, rwdlow));
    set(gca,'ydir','reverse', 'xdir','reverse')  

% Proportion savings barplots for effort scale    
nexttile(3, [3 1])
    b = barh(categorical({'Young' 'Old'}), [effpropsyoung(2,3)' effpropsyoung(1,3)'; effpropsold(2,3)' effpropsold(1,3)'],'stacked');
    b(1).FaceColor = 'flat';
    b(1).CData = [2/5 0 0; 0 2/5 0;];
    b(2).FaceColor = 'flat';
    b(2).CData = [5/5 0 0; 0 5/5 0];
    ylabel(sprintf('Increased effort'));
    xlabel(sprintf(['Proportion of savings\n%dJ - %dJ'], rwdhigh, rwdlow));

% Proportion savings barplots for reward scale    
nexttile(12, [3 1])
    b = barh(categorical({'High' 'Low'}), [rwdpropsyoung(2,3)' rwdpropsyoung(1,3)'; rwdpropsyoung(2,1)' rwdpropsyoung(1,1)'],'stacked');
    b(1).FaceColor = 'flat';
    b(1).CData = [0 2/5 0; 2/5 0 0];
    b(2).FaceColor = 'flat';
    b(2).CData = [0 5/5 0; 5/5 0 0];
    ylabel(sprintf('Reward scale'));
    xlabel(sprintf(['Proportion of savings\n%dJ - %dJ'], rwdhigh, rwdlow));
    set(gca,'ydir','reverse')

% Proportion savings barplots for probability scale
nexttile(21, [3 1])
    b = barh(categorical({'High' 'Low'}), [probpropsyoung(2,3)' probpropsyoung(1,3)'; probpropsyoung(2,5)' probpropsyoung(1,5)'],'stacked');
    b(1).FaceColor = 'flat';
    b(1).CData = [0 2/5 0; 2/5 0 0];
    b(2).FaceColor = 'flat';
    b(2).CData = [0 5/5 0; 5/5 0 0];
    ylabel(sprintf('Reward scale'));
    xlabel(sprintf(['Proportion of savings\n%dJ - %dJ'], rwdhigh, rwdlow));
    set(gca,'ydir','reverse')

beautifyfig;
