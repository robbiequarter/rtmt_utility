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
rwdlow = 50;
rwdhigh = 55;
alphalow = find(myalphas==rwdlow);
alphahigh = find(myalphas==rwdhigh);

%% Process data

% Organize model data
for j=1:length(myeffscales)
    
    effcurvesyoung=squeeze(young(:,j,alphascaleind,probscaleind,:));
    effcurvesyoung_mt(:,j)=effcurvesyoung(:,1);
    effcurvesyoung_rt(:,j)=effcurvesyoung(:,2);
    effsavings = ((effcurvesyoung_rt(alphahigh,j) - effcurvesyoung_rt(alphalow,j))+(effcurvesyoung_mt(alphahigh,j) - effcurvesyoung_mt(alphalow,j))); % total time savings
    effdeltasyoung(1,j) = (effcurvesyoung_mt(alphahigh,j) - effcurvesyoung_mt(alphalow,j));
    effdeltasyoung(2,j) = (effcurvesyoung_rt(alphahigh,j) - effcurvesyoung_rt(alphalow,j));
    effpropsyoung(1,j) = (effcurvesyoung_mt(alphahigh,j) - effcurvesyoung_mt(alphalow,j))/effsavings; % MT proportion of time savings
    effpropsyoung(2,j) = (effcurvesyoung_rt(alphahigh,j) - effcurvesyoung_rt(alphalow,j))/effsavings; % RT proportion of time savings
    stackData1(j,2,1) = effpropsyoung(1,j);
    stackData1(j,2,2) = effpropsyoung(2,j);
end

for k=1:length(myalphascales)
    
    alphacurvesyoung=squeeze(young(:,effscaleind,k,probscaleind,:));
    alphacurvesyoung_mt(:,k)=alphacurvesyoung(:,1);
    alphacurvesyoung_rt(:,k)=alphacurvesyoung(:,2);
    rwdsavings = ((alphacurvesyoung_rt(alphahigh,k) - alphacurvesyoung_rt(alphalow,k))+(alphacurvesyoung_mt(alphahigh,k) - alphacurvesyoung_mt(alphalow,k)));
    rwddeltasyoung(1,k) = (alphacurvesyoung_mt(alphahigh,k) - alphacurvesyoung_mt(alphalow,k));
    rwddeltasyoung(2,k) = (alphacurvesyoung_rt(alphahigh,k) - alphacurvesyoung_rt(alphalow,k));
    rwdpropsyoung(1,k) = (alphacurvesyoung_mt(alphahigh,k) - alphacurvesyoung_mt(alphalow,k))/rwdsavings; % MT proportion of time savings
    rwdpropsyoung(2,k) = (alphacurvesyoung_rt(alphahigh,k) - alphacurvesyoung_rt(alphalow,k))/rwdsavings; % RT proportion of time savings
    stackData2(k,2,1) = rwdpropsyoung(1,k);
    stackData2(k,2,2) = rwdpropsyoung(2,k);
end

% Proportions of savings from empirical data
dattimes(:,1) = table2array(aggdat(aggdat.MASS ==0 & aggdat.REWD == 1, [14,16]))' - table2array(aggdat(aggdat.MASS ==0 & aggdat.REWD == 0, [14,16]))'; % col 1 = no mass, row 1 = RT
dattimes(:,2) = table2array(aggdat(aggdat.MASS ==1 & aggdat.REWD == 1, [14,16]))' - table2array(aggdat(aggdat.MASS ==1 & aggdat.REWD == 0, [14,16]))';
datprops(1,1) = dattimes(2,1)/sum(dattimes(:,1)); datprops(2,1) = 1-datprops(1,1);
datprops(1,2) = dattimes(2,2)/sum(dattimes(:,2)); datprops(2,2) = 1-datprops(1,2);


%% Proprotion barplots, separated MT/RT
figure
    subplot(2,2,1)
        bar(categorical(myeffscales), effpropsyoung(2,:)','r');
    %     ylim([-0.07 0]);
        xlabel('Increasing Effort Scaling');
        ylabel(sprintf(['Proportion of savings\n%dJ - %dJ'], rwdhigh, rwdlow));
        title('Reaction Time');
    subplot(2,2,2)
        bar(categorical(myeffscales), effpropsyoung(1,:)');
    %     ylim([-0.07 0]);
        xlabel('Increasing Effort Scaling');
        ylabel(sprintf(['Proportion of savings\n%dJ - %dJ'], rwdhigh, rwdlow));
        title('Movement Time')
     subplot(2,2,3)
        bar(categorical(myalphascales), rwdpropsyoung(2,:),'r');
        set(gca,'xdir','reverse')
    %     ylim([-0.07 0]);
        xlabel('Reducing Reward Scaling');
        ylabel(sprintf(['Proportion of savings\n%dJ - %dJ'], rwdhigh, rwdlow));
        title('Reaction Time');
    subplot(2,2,4)
        bar(categorical(myalphascales), rwdpropsyoung(1,:));
    %     ylim([-0.07 0]);
        set(gca,'xdir','reverse')
        xlabel('Reducing Reward Scaling');
        ylabel(sprintf(['Proportion of savings\n%dJ - %dJ'], rwdhigh, rwdlow));
        title('Movement Time')

%% Figure for horizontal proportions, stacked
figure
    subplot(4,1,1)
        [tme, r, v, a] = minjerk([0,0],[0,0.1],1,.001);
        plot(tme,v(:,2))
    subplot(4,1,2)
        barh(categorical(myeffscales), [effpropsyoung(2,:)' effpropsyoung(1,:)'],'stacked');
        ylabel(sprintf('Increasing Effort\nScale'));
        xlabel(sprintf(['Proportion of savings\n%dJ - %dJ'], rwdhigh, rwdlow));
        set(gca,'ydir','reverse');
        legend('RT','MT','Location','southeast');
    subplot(4,1,3)
        barh(categorical(myalphascales), [rwdpropsyoung(2,:)' rwdpropsyoung(1,:)'] ,'stacked');
        ylabel(sprintf('Reducing Reward\nScale'));
        xlabel(sprintf(['Proportion of savings\n%dJ - %dJ'], rwdhigh, rwdlow));
    subplot(4,1,4)
        barh(categorical({'LOW','HIGH'}), [datprops(2,:)' datprops(1,:)'] ,'stacked');
        ylabel('Effort');
        xlabel(sprintf(['Proportion of savings'], rwdhigh, rwdlow));
    beautifyfig;
    
%% Re-worked curves, x-axis is scaling factor, each line is different rewards
% alphas = 21:10:51;
% for j = 1:length(myeffscales)
%         effmt(:,j) = squeeze(young(alphas,j,alphascaleind,probscaleind,1)); % MT 
%         effrt(:,j) = squeeze(young(alphas,j,alphascaleind,probscaleind,2)); 
% 
%         rwdmt(:,j) = squeeze(young(alphas,effscaleind,j,probscaleind,1)); 
%         rwdrt(:,j) = squeeze(young(alphas,effscaleind,j,probscaleind,2));
% end
% 
% 
% 
% figure
% for j = 1:length(alphas)
%         subplot(2,1,1)          
%             h(1,1)=plot(myeffscales,effmt(j,:),'o-'); % MT
%             hold on 
%             h(1,2)=plot(myeffscales,effrt(j,:),'x-'); % RT
%             %h(1,2).LineStyle = '--'; % RT = dashed
%             % Thicker and more colored = Increasing effort scale
%                 h(1,1).Color=j*[1 0 0]./length(alphas);
%                 h(1,2).Color=j*[0 0 1]./length(alphas);
%                 h(1,1).LineWidth = 1; 
%                 h(1,2).LineWidth = 1;
%             legend(h(1,:),'MT','RT')
%             xlabel('Effort Scale')
%             ylabel('Duration (s)')
% %             set(gca,'ylim',[0.3 1.4])
% 
%         subplot(2,1,2)          
%             h(1,1)=plot(myeffscales,rwdmt(j,:),'o-'); % MT
%             hold on 
%             h(1,2)=plot(myeffscales,rwdrt(j,:),'x-'); % RT
%             %h(1,2).LineStyle = '--'; % RT = dashed
%             % Thicker and more colored = Increasing effort scale
%                 h(1,1).Color=j*[1 0 0]./length(alphas);
%                 h(1,2).Color=j*[0 0 1]./length(alphas);
%                 h(1,1).LineWidth = 1; 
%                 h(1,2).LineWidth = 1;
%             legend('MT','RT')
%             xlabel('Reward Scale')
%             ylabel('Duration (s)')
%             set(gca,'xdir','reverse');
% end









