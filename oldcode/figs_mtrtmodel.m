%% Re-creating figures with simulated old/young data
% clear all
% close all

%% 1. Load in data for each simulation
% simresults_old2/young2 is only scaling "b" term
% simresults_old3 only scales up old "a" by 1.6
% simresults_old4 changes old resting rate (ao) from 69 to 77 to match young
% simresults_old5/young5 is new RT function (c1 = 35; c0 = -8.75)
% simresults_old6/young6 is new RT function (c1 = 20; c0 = -6)

% young
load('simresults/simresults_young.mat'); 
young = mysols;

% old
load('simresults/simresults_old.mat'); 
old = mysols;

%% 2. Set parameters
global param

%metabolic data parameters: etotal = esit + emove = ao*to + a*t + (b*d^i)/T

%old
param.myc0 = -9; % accuracy parameters; shifts logistic to the right with scaling
param.myc1 = 15; % accuracy parameters; shifts logistic to the left with scaling
param.mya = 77; % effort offset, try scaling this up by 1.6
param.myb = 20; %b in metabolic equation 
param.myao = 69; % resting rate for older adults, try matching this to young also
% param.myao = 77; % resting rate for older adults, try matching this to young also
param.myi=0.88;  %exponent on distance

%young
param.myc0young = -5;  
param.myc1young = 10; 
param.myayoung = 77;
param.mybyoung = 12;
param.myaoyoung = 77;
param.myiyoung=1.23;

%both
param.myc0to = -1; % Used in optimization
param.myc1to = 10; % Used in optimization
% param.myc0to = -6; % -8.75; % Used in new RT function
% param.myc1to = 20; % 35; % Used in new RT function
param.myeffscale = 1; 

%distance
d=0.1;

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

%% 3. Effort Valuation Figure #1
%Effort Valuation Scaling

% figure
% for j = 1:length(myeffscales)
%     subplot(1,5,j)    
%         h(:,1)=plot(myalphas,squeeze(old(:,j,alphascaleind,probscaleind,:)),...
%             'Color', 'r', 'LineWidth', 1);
%         hold on
%         h(2,1).LineStyle = '--';
%         h(:,2)=plot(myalphas,squeeze(young(:,j,alphascaleind,probscaleind,:)),...
%             'Color', 'g', 'LineWidth',1);
%         h(2,2).LineStyle = '--';
% 
%         title(sprintf(['Increasing effort valuation\n (effort scaling)']))
%         legend(h(1,:),'Old','Young')
%         xlabel('Reward (alpha,J)')
%         ylabel('Duration (s)')
%         set(gca,'ylim',[0.3 1.4])
%         hold off
% end

%% 4. Effort & Reward Valuation Figure
figure
for j = 1:length(myeffscales)
    subplot(2,2,1)    
        h(:,1)=plot(myalphas,squeeze(old(:,j,alphascaleind,probscaleind,:)));
        hold on
        h(2,1).LineStyle = '--'; % RT = dashed
        % Thicker and more colored = Increasing effort scale
            h(1,1).Color=j*[1 0 0]./length(myeffscales); 
            h(2,1).Color=j*[1 0 0]./length(myeffscales);
            h(1,1).LineWidth=1; %0.3*j
            h(2,1).LineWidth=1;
        legend(h(:,1),'MT','RT')
        xlabel('Reward (alpha,J)')
        ylabel('Duration (s)')
        title('OLD')
        set(gca,'ylim',[0.3 1.4])
        
    subplot(2,2,2)
        h(:,2)=plot(myalphas,squeeze(young(:,j,alphascaleind,probscaleind,:)));
        hold on 
        h(2,2).LineStyle = '--'; % RT = dashed
        % Thicker and more colored = Increasing effort scale
            h(1,2).Color=j*[0 1 0]./length(myeffscales);
            h(2,2).Color=j*[0 1 0]./length(myeffscales);
            h(1,2).LineWidth = 1; 
            h(2,2).LineWidth = 1;
        legend(h(:,2),'MT','RT')
        xlabel('Reward (alpha,J)')
        ylabel('Duration (s)')
        title('YOUNG')
        set(gca,'ylim',[0.3 1.4])
end
% sgtitle(sprintf(['Increasing effort valuation\n (effscaling)']));

%Reward Valuation Scaling
for k = 1:length(myalphascales)
    subplot(2,2,3)    
        h(:,1)=plot(myalphas,squeeze(old(:,effscaleind,k,probscaleind,:)));
        hold on
        h(2,1).LineStyle = '--'; % RT = dashed
        % Thicker and more colored = Increasing reward scale
            h(1,1).Color=k*[1 0 0]./length(myalphascales); 
            h(2,1).Color=k*[1 0 0]./length(myalphascales);
            h(1,1).LineWidth=1;
            h(2,1).LineWidth=1;
        legend(h(:,1),'MT','RT')
        xlabel('Reward (alpha,J)')
        ylabel('Duration (s)')
        title('OLD')
        set(gca,'ylim',[0.3 1.4])
        
    subplot(2,2,4)
        h(:,2)=plot(myalphas,squeeze(young(:,effscaleind,k,probscaleind,:)));
        hold on 
        h(2,2).LineStyle = '--'; % RT = dashed
        % Thicker and more colored = Increasing effort scale
            h(1,2).Color=k*[0 1 0]./length(myalphascales);
            h(2,2).Color=k*[0 1 0]./length(myalphascales);
            h(1,2).LineWidth=1; 
            h(2,2).LineWidth=1;
        legend(h(:,2),'MT','RT')
        xlabel('Reward (alpha,J)')
        ylabel('Duration (s)')
        title('YOUNG')
        set(gca,'ylim',[0.3 1.4])
end
sgtitle(sprintf(['Top Row: Effort Scaling\n Bottom Row: Alpha Scaling']));

%% 5. Re-create the Hackier Figure
figure
j = 3; %effindex
h(:,1)=plot(myalphas,squeeze(old(:,j,alphascaleind,probscaleind,:)),...
            'Color', 'r', 'LineWidth', 1);
hold on
h(2,1).LineStyle = '--';
h(:,2)=plot(myalphas,squeeze(young(:,j,alphascaleind,probscaleind,:)),...
    'Color', 'g', 'LineWidth',1);
h(2,2).LineStyle = '--';

title(sprintf(['Increasing effort valuation\n (effort scale = %.1f)'],myeffscales(j)))
legend(h(1,:),'Old','Young')
xlabel('Reward (alpha,J)')
ylabel('Duration (s)')
set(gca,'ylim',[0.3 1.4])
beautifyfig;

%% 6. Delta RT/MT for a Change in Reward Barplots
% Maybe try like reward = 100 - 50 for different effort scalings
rwdlow = 50;
rwdhigh = 55;
alphalow = find(myalphas==rwdlow);
alphahigh = find(myalphas==rwdhigh);

for j=1:length(myeffscales)
    effcurvesold=squeeze(old(:,j,alphascaleind,probscaleind,:));
    effcurvesold_mt(:,j)=effcurvesold(:,1);
    effcurvesold_rt(:,j)=effcurvesold(:,2);   
    effdeltasold(1,j) = effcurvesold_mt(alphahigh,j) - effcurvesold_mt(alphalow,j);
    effdeltasold(2,j) = effcurvesold_rt(alphahigh,j) - effcurvesold_rt(alphalow,j);
    stackData1(j,1,1) = effdeltasold(1,j);
    stackData1(j,1,2) = effdeltasold(2,j);
    
    effcurvesyoung=squeeze(young(:,j,alphascaleind,probscaleind,:));
    effcurvesyoung_mt(:,j)=effcurvesyoung(:,1);
    effcurvesyoung_rt(:,j)=effcurvesyoung(:,2);
    effdeltasyoung(1,j) = effcurvesyoung_mt(alphahigh,j) - effcurvesyoung_mt(alphalow,j);
    effdeltasyoung(2,j) = effcurvesyoung_rt(alphahigh,j) - effcurvesyoung_rt(alphalow,j);
    stackData1(j,2,1) = effdeltasyoung(1,j);
    stackData1(j,2,2) = effdeltasyoung(2,j);
end

for k=1:length(myalphascales)
    alphacurvesold=squeeze(old(:,effscaleind,k,probscaleind,:));
    alphacurvesold_mt(:,k)=alphacurvesold(:,1);
    alphacurvesold_rt(:,k)=alphacurvesold(:,2);
    rwddeltasold(1,k) = alphacurvesold_mt(alphahigh,k) - alphacurvesold_mt(alphalow,k);
    rwddeltasold(2,k) = alphacurvesold_rt(alphahigh,k) - alphacurvesold_rt(alphalow,k);
    stackData2(k,1,1) = rwddeltasold(1,k);
    stackData2(k,1,2) = rwddeltasold(2,k);
    
    alphacurvesyoung=squeeze(young(:,effscaleind,k,probscaleind,:));
    alphacurvesyoung_mt(:,k)=alphacurvesyoung(:,1);
    alphacurvesyoung_rt(:,k)=alphacurvesyoung(:,2);
    rwddeltasyoung(1,k) = alphacurvesyoung_mt(alphahigh,k) - alphacurvesyoung_mt(alphalow,k);
    rwddeltasyoung(2,k) = alphacurvesyoung_rt(alphahigh,k) - alphacurvesyoung_rt(alphalow,k);
    stackData2(k,2,1) = rwddeltasyoung(1,k);
    stackData2(k,2,2) = rwddeltasyoung(2,k);
end

figure
subplot(2,2,1)
    bar(myeffscales, abs(effdeltasold'));
%     ylim([0 0.8]);
    xlabel('Effort Scaling');
    ylabel(sprintf(['\\DeltaT (s)\n%dJ - %dJ'], rwdhigh, rwdlow));
    title('Old');
    legend('MT', 'RT','Location', 'southeast');
subplot(2,2,2)
    bar(myeffscales, abs(effdeltasyoung'));
%     ylim([0 0.8]);
    xlabel('Effort Scaling');
    ylabel(sprintf(['\\DeltaT (s)\n%dJ - %dJ'], rwdhigh, rwdlow));
    title('Young')
 subplot(2,2,3)
    bar(myalphascales, abs(rwddeltasold'));
%     ylim([0 0.8]);
    xlabel('Reward Scaling');
    ylabel(sprintf(['\\DeltaT (s)\n%dJ - %dJ'], rwdhigh, rwdlow));
    title('Old');
    legend('MT', 'RT','Location', 'southeast');
subplot(2,2,4)
    bar(myalphascales, abs(rwddeltasyoung'));
%     ylim([0 0.8]);
    xlabel('Reward Scaling');
    ylabel(sprintf(['\\DeltaT (s)\n%dJ - %dJ'], rwdhigh, rwdlow));
    title('Young')
    

%% Scatter plots of deltas with best fit line 
figure
subplot(2,1,1)
    % Effort scaling
    plot(myeffscales, abs(effdeltasold), 'r.','MarkerSize', 10);
    hold on
    plot(myeffscales, abs(effdeltasyoung), 'g.','MarkerSize', 10);
    h2 = lsline; % for some reason, lsline puts RT first
    % I think lsline might work backwards from the order of what was
    % plotted
    h2(1,1).LineStyle = '--'; h2(1,3).LineStyle = '--';
    hold off 
    ylim([0 0.1]);
%     legend(h2,'Young \Delta RT','Young \Delta MT','Old \Delta RT','Old \Delta MT')
    xlabel('Effort Scaling');
    ylabel(sprintf(['\\DeltaT (s)\n%dJ - %dJ'], rwdhigh, rwdlow));
    
subplot(2,1,2)
    % RWD scaling
    plot(myeffscales, abs(rwddeltasold), 'r.','MarkerSize', 10);
    hold on
    plot(myeffscales, abs(rwddeltasyoung), 'g.','MarkerSize', 10);
    h2 = lsline; % for some reason, lsline puts RT first
    % I think lsline might work backwards from the order of what was
    % plotted
    h2(1,1).LineStyle = '--'; h2(1,3).LineStyle = '--';
    hold off 
    ylim([0 0.1]);
%     legend(h2,'Young \Delta RT','Young \Delta MT','Old \Delta RT','Old \Delta MT')
    xlabel('Reward Scaling');
    ylabel(sprintf(['\\DeltaT (s)\n%dJ - %dJ'], rwdhigh, rwdlow));
    
%% Scatter Plots w/ Sliding Reward Diffs
rwdlows = 50:1:70;
rwdhighs = 51:1:71;

count = 1;
for i = 1:length(rwdlows)
    alphalow = find(myalphas==rwdlows(i));
    alphahigh = find(myalphas==rwdhighs(i));
    for j=1:length(myeffscales)
        effcurvesold=squeeze(old(:,j,alphascaleind,probscaleind,:));
        effcurvesold_mt(:,j)=effcurvesold(:,1);
        effcurvesold_rt(:,j)=effcurvesold(:,2);   
        effdeltas_oldmt(count,1) = myeffscales(j);
        effdeltas_oldmt(count,2) = effcurvesold_mt(alphahigh,j) - effcurvesold_mt(alphalow,j);
        effdeltas_oldrt(count,1) = myeffscales(j);
        effdeltas_oldrt(count,2) = effcurvesold_rt(alphahigh,j)  - effcurvesold_rt(alphalow,j);
        
        effcurvesyoung=squeeze(young(:,j,alphascaleind,probscaleind,:));
        effcurvesyoung_mt(:,j)=effcurvesyoung(:,1);
        effcurvesyoung_rt(:,j)=effcurvesyoung(:,2);
        effdeltas_youngmt(count,1) = myeffscales(j);
        effdeltas_youngmt(count,2) = effcurvesyoung_mt(alphahigh,j) - effcurvesyoung_mt(alphalow,j);
        effdeltas_youngrt(count,1) = myeffscales(j);
        effdeltas_youngrt(count,2) = effcurvesyoung_rt(alphahigh,j)  - effcurvesyoung_rt(alphalow,j);
        
        alphacurvesold=squeeze(old(:,effscaleind,k,probscaleind,:));
        alphacurvesold_mt(:,j)=alphacurvesold(:,1);
        alphacurvesold_rt(:,j)=alphacurvesold(:,2);        
        rwddeltas_oldmt(count,1) = myeffscales(j);
        rwddeltas_oldmt(count,2) = alphacurvesold_mt(alphahigh,j) - alphacurvesold_mt(alphalow,j);
        rwddeltas_oldrt(count,1) = myeffscales(j);
        rwddeltas_oldrt(count,2) = alphacurvesold_rt(alphahigh,j)  - alphacurvesold_rt(alphalow,j);
        
        alphacurvesyoung=squeeze(young(:,effscaleind,k,probscaleind,:));
        alphacurvesyoung_mt(:,j)=alphacurvesyoung(:,1);
        alphacurvesyoung_rt(:,j)=alphacurvesyoung(:,2);
        rwddeltas_youngmt(count,1) = myeffscales(j);
        rwddeltas_youngmt(count,2) = alphacurvesyoung_mt(alphahigh,j) - alphacurvesyoung_mt(alphalow,j);
        rwddeltas_youngrt(count,1) = myeffscales(j);
        rwddeltas_youngrt(count,2) = alphacurvesyoung_rt(alphahigh,j)  - alphacurvesyoung_rt(alphalow,j);
        count = count+1;
    end
end

efffit_mtold = polyfit(effdeltas_oldmt(:,1),effdeltas_oldmt(:,2),1);
efffit_rtold = polyfit(effdeltas_oldrt(:,1),effdeltas_oldrt(:,2),1);
efffit_mtyoung = polyfit(effdeltas_youngmt(:,1),effdeltas_youngmt(:,2),1);
efffit_rtyoung = polyfit(effdeltas_youngrt(:,1),effdeltas_youngrt(:,2),1);

rwdfit_mtold = polyfit(rwddeltas_oldmt(:,1),rwddeltas_oldmt(:,2),1);
rwdfit_rtold = polyfit(rwddeltas_oldrt(:,1),rwddeltas_oldrt(:,2),1);
rwdfit_mtyoung = polyfit(rwddeltas_youngmt(:,1),rwddeltas_youngmt(:,2),1);
rwdfit_rtyoung = polyfit(rwddeltas_youngrt(:,1),rwddeltas_youngrt(:,2),1);

figure
subplot(1,2,1)
    % Effort scaling - old
    plot(effdeltas_oldmt(:,1), effdeltas_oldmt(:,2), 'r.', 'MarkerSize', 10);
    hold on
    plot(effdeltas_oldrt(:,1), effdeltas_oldrt(:,2), 'r.','MarkerSize', 10);
    plot(myeffscales, polyval(efffit_mtold,myeffscales), 'r', 'LineWidth',1)
    plot(myeffscales, polyval(efffit_rtold,myeffscales), 'r--', 'LineWidth',1)

    % Effort scaling - young
    plot(effdeltas_youngmt(:,1), effdeltas_youngmt(:,2), 'g.','MarkerSize', 10);
    plot(effdeltas_youngrt(:,1), effdeltas_youngrt(:,2), 'g.','MarkerSize', 10);
    plot(myeffscales, polyval(efffit_mtyoung,myeffscales), 'g', 'LineWidth',1)
    plot(myeffscales, polyval(efffit_rtyoung,myeffscales), 'g--', 'LineWidth',1)
    hold off
    xlabel('Effort Scaling');
    ylabel('\DeltaT (s)');
    
subplot(1,2,2)
    % RWD scaling - old
    plot(rwddeltas_oldmt(:,1), rwddeltas_oldmt(:,2), 'r.', 'MarkerSize', 10);
    hold on
    plot(rwddeltas_oldrt(:,1), rwddeltas_oldrt(:,2), 'r.','MarkerSize', 10);
    plot(myalphascales, polyval(rwdfit_mtold,myalphascales), 'r', 'LineWidth',1)
    plot(myalphascales, polyval(rwdfit_rtold,myalphascales), 'r--', 'LineWidth',1)

    % RWD scaling - young
    plot(rwddeltas_youngmt(:,1), rwddeltas_youngmt(:,2), 'g.','MarkerSize', 10);
    plot(rwddeltas_youngrt(:,1), rwddeltas_youngrt(:,2), 'g.','MarkerSize', 10);
    plot(myalphascales, polyval(rwdfit_mtyoung,myalphascales), 'g', 'LineWidth',1)
    plot(myalphascales, polyval(rwdfit_rtyoung,myalphascales), 'g--', 'LineWidth',1)
    hold off
    
    xlabel('Reward Scaling');
    ylabel('\DeltaT (s)');

%% 7. Grouped + Stacked Barplots
% Uses plotBarStackGroups.m from File Exchange:
% Evan (2021). Plot Groups of Stacked Bars (https://www.mathworks.com/matlabcentral/fileexchange/32884-plot-groups-of-stacked-bars), 
%   MATLAB Central File Exchange. Retrieved May 18, 2021. 

figure
groupLabels = {'0.8', '0.9', '1.0', '1.1', '1.2'};
subplot(2,1,1)
    plotBarStackGroups(stackData1, groupLabels)
    xlabel('Effort Scaling');
    ylabel(sprintf(['\\DeltaT (s)\n%dJ - %dJ'], rwdhigh, rwdlow));
%     ylim([-0.8, 0]);
    legend('MT Old', 'RT Old', 'MT Young', 'RT Young',...
        'Location', 'south', 'Orientation', 'horizontal');
subplot(2,1,2)
    plotBarStackGroups(stackData2, groupLabels)
    xlabel('Reward Scaling');
    ylabel(sprintf(['\\DeltaT (s)\n%dJ - %dJ'], rwdhigh, rwdlow));
%     ylim([-0.8, 0]);
    legend('MT Old', 'RT Old', 'MT Young', 'RT Young',...
        'Location', 'south', 'Orientation', 'horizontal');
    
%% 8. Playing around with RT Probability of Reward function
% Logisitic curve: f(x) = L/(1+exp(-k(x-x0))) where:
%   L = function max, so 1 for probability
%   x0 = sigmoid (halfway) point. In new function, this is at 0.625 prob. 
%   k = growth rate (steepness)
% Ideally want anything under 150ms to be chance, so 25% hitting the target
    % Might actually be around 100-150ms (Haith et al., 2016)
    % Looks like the 100ms (Haith et al., 2016) was in forced RT, 150ms was
    % free RT
% Solution: ln(3) = k*x0 for P(0) = 0.25
% Rapid increase to 1 afterwards
param.myc0to = -1;
param.myc1to = 10;

% k = 30, x0 = 0.25 
% k = 25, x0 = 0.35
% k = 20, x0 = 0.45
% k = 20; % same as c1
% x0 = 0.45; % -1*c0/c1
k = 35; % same as c1
x0 = 0.25; % -1*c0/c1
k00 = 20; 
x00 = 0.30;

Pr_test = @(Tr) (0.75./(1 + exp(-k.*(Tr - x0))))+0.25;
Pr_test2 = @(Tr) (0.75./(1 + exp(-k00.*(Tr - x00))))+0.25;
Pr_oldyoung = @(Tr) (1./(1 + exp(-(param.myc0to) - (param.myc1to).*Tr)));
% Pr_oldyoung2 = @(Tr) (0.75./(1 + exp(-(x00) - (k00).*Tr)))+0.25;

figure
    fplot(Pr_oldyoung,'r', 'LineWidth', 1, 'DisplayName', 'old rt function')
    hold on
    % fplot(Pr_oldyoung2, 'r--', 'LineWidth', 1)
    fplot(Pr_test, 'Color', [0 2/6 0], 'LineWidth', 1, 'DisplayName',sprintf('k = %d; x_0 = %.2f',k,x0))
    fplot(Pr_test2, 'Color', [0 4/6 0], 'LineWidth', 1, 'DisplayName',sprintf('k = %d; x_0 = %.2f',k00,x00))
    xline(0.150,'k--','HandleVisibility','off'); yline(.25,'k--','HandleVisibility','off');
    % xline(x0,'k--'); yline(.5,'k--');
    legend('show', 'Location','southeast');
    xlim([-0.5 1.0]);
    % title(sprintf(['k = %d; x_0 = %.2f'],k,x0));
    hold off

k2 = 35:-5:10;
x02 = 0.15:0.05:0.4;
% x02 = 0.25.*ones(size(k2));
% k2 = 25.*ones(size(x02));
figure
legend_names = cell(1,length(k2));
for j = 1:length(k2)
    k = k2(j);
    x0 = x02(j);
    Pr_test = @(Tr) (0.75./(1 + exp(-k.*(Tr - x0))))+0.25;
    z(1,1) = fplot(Pr_test,'LineWidth',1);
    hold on
    z(1,1).Color=j*[0 1 0]./length(x02);
    legend_names{1,j} = sprintf(['k = %d; x_0 = %.2f'],k,x0);
end
fplot(Pr_oldyoung,'r', 'LineWidth', 1)
xline(0.150,'k--', 'LineWidth',1); yline(.25,'k--','LineWidth',1);
xlim([-0.5 1.0]); ylim([0, 1]);
title(sprintf('More Green = Smaller k, Larger x_0\nPr = 0.75/(1+exp(-k(Tr - x0))) + 0.25'));
xlabel('RT (s)'); ylabel('Pr(RWD)');
legend(legend_names,'Location','southeast');
caption = sprintf('Pr = 0.75/(1+exp(-k(Tr-x0))) + 0.25');
hold off