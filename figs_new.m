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
load('simresults/simresults_young10.mat'); 
young = mysols;
% old
load('simresults/simresults_old10.mat'); 
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

%% Proprotion barplots, separated MT/RT
% figure
%     subplot(2,2,1)
%         bar(categorical(myeffscales), effpropsyoung(2,:)','r');
%     %     ylim([-0.07 0]);
%         xlabel('Increasing Effort Scaling');
%         ylabel(sprintf(['Proportion of savings\n%dJ - %dJ'], rwdhigh, rwdlow));
%         title('Reaction Time');
%     subplot(2,2,2)
%         bar(categorical(myeffscales), effpropsyoung(1,:)');
%     %     ylim([-0.07 0]);
%         xlabel('Increasing Effort Scaling');
%         ylabel(sprintf(['Proportion of savings\n%dJ - %dJ'], rwdhigh, rwdlow));
%         title('Movement Time')
%      subplot(2,2,3)
%         bar(categorical(myalphascales), rwdpropsyoung(2,:),'r');
%         set(gca,'xdir','reverse')
%     %     ylim([-0.07 0]);
%         xlabel('Reducing Reward Scaling');
%         ylabel(sprintf(['Proportion of savings\n%dJ - %dJ'], rwdhigh, rwdlow));
%         title('Reaction Time');
%     subplot(2,2,4)
%         bar(categorical(myalphascales), rwdpropsyoung(1,:));
%     %     ylim([-0.07 0]);
%         set(gca,'xdir','reverse')
%         xlabel('Reducing Reward Scaling');
%         ylabel(sprintf(['Proportion of savings\n%dJ - %dJ'], rwdhigh, rwdlow));
%         title('Movement Time')

%% Figure for horizontal proportions, stacked
figure
%     subplot(4,1,1)
%         [tme, r, v, a] = minjerk([0,0],[0,0.1],1,.001);
%         plot(tme,v(:,2))
%     subplot(4,1,2)
    subplot(3,1,1)
        barh(categorical(myeffscales), [effpropsyoung(2,:)' effpropsyoung(1,:)'],'stacked');
        ylabel(sprintf('Increasing Effort\nScale'));
        xlabel(sprintf(['Proportion of savings\n%dJ - %dJ'], rwdhigh, rwdlow));
        set(gca,'ydir','reverse');
        legend('RT','MT','Location','southeast');
%     subplot(4,1,3)
    subplot(3,1,2)
        barh(categorical(myalphascales), [rwdpropsyoung(2,:)' rwdpropsyoung(1,:)'] ,'stacked');
        ylabel(sprintf('Reducing Reward\nScale'));
        xlabel(sprintf(['Proportion of savings\n%dJ - %dJ'], rwdhigh, rwdlow));
%     subplot(4,1,4)
    subplot(3,1,3)
        x = categorical({'LOW','HIGH' 'YOUNG' 'OLD'});
        barh(reordercats(x, string(x)), [datprops(2,:)' datprops(1,:)'] ,'stacked');
        ylabel('Effort');
        xlabel(sprintf('Proportion of savings\nRWD - NRWD (a.u.)'));
    beautifyfig;

%% Combined figure

figure
clear h
% ADD 50-55J SHADING BOX
% 2x2 for first figure, 2x1 for RT/MT bars, 2x2 for last line graphs

% for shaded/highlighted patch of plot
ptchidx = (myalphas >= rwdlow) & (myalphas <= rwdhigh);
ptchcol = [221,160,221]./256;
colgray = [0.6    0.6    0.6];
unit = 1000; %scaling factor (to go from s to ms)

% subplot(4,5,[1,2,6,7])
subplot(4,2,[1,3])
    for i = 1:5
        h(:,1)=plot(myalphas,squeeze(young(:,i,alphascaleind,probscaleind,:)).*unit,'LineWidth', 1);
        hold on
        h(1,1).Color = 1/i.*[0 1 0];
        h(2,1).Color = 1/i.*[0 1 0];
        h(2,1).LineStyle = '--';
    %     h(2,1).Color = 'r';
%         h(:,2)=plot(myalphas,squeeze(young(:,5,alphascaleind,probscaleind,:)).*unit,...
%             'Color', 'r', 'LineWidth',1);
%         h(2,2).LineStyle = '--';
    %     h(2,2).Color = 'r';
    %     title('Increasing effort valuation')
    end
    patch([myalphas(ptchidx) fliplr(myalphas(ptchidx))], [1.5.*ones(size(myalphas(ptchidx))).*unit -0.5.*ones(size(myalphas(ptchidx))).*unit],...
    ptchcol, 'FaceAlpha',0.3, 'EdgeColor','none')
    legend('EffScale = 0.8')
%     xlabel('Reward (alpha,J)')
    xticks([40 60 80 100])
    yticks([300 600 900 1200 1500])
    ylabel('Duration (ms)')
    set(gca,'ylim',[0.3 1.2].*unit,'xlim',[35 65])
    
% subplot(4,5,[11,12,16,17])
subplot(4,2,[5,7])
    for i = 1:5
        h(:,1)=plot(myalphas,squeeze(young(:,effscaleind,i,probscaleind,:)).*unit,...
                    'LineWidth', 1);
        hold on
        h(2,1).LineStyle = '--';
        h(1,1).Color = i.*[1/5 0 0];
        h(2,1).Color = i.*[1/5 0 0];
    end
    patch([myalphas(ptchidx) fliplr(myalphas(ptchidx))], [1.5.*ones(size(myalphas(ptchidx))).*unit -0.5.*ones(size(myalphas(ptchidx))).*unit],...
        ptchcol, 'FaceAlpha',0.3, 'EdgeColor','none')
    legend('RWDScale = 0.8')
    xlabel('Reward (alpha,J)')
    ylabel('Duration (ms)')
    xticks([40 60 80 100])
    yticks([300 600 900 1200 1500])
    set(gca,'ylim',[0.3 1.2].*unit,'xlim',[35 65])
    
% subplot(4,5,3)
subplot(4,2,[2 4])
    barh(categorical(myeffscales), [effpropsyoung(2,:)' effpropsyoung(1,:)'],'stacked');
    ylabel(sprintf('Increasing Effort\nScale'));
%     xlabel(sprintf(['Proportion of savings\n%dJ - %dJ'], rwdhigh, rwdlow));
    set(gca,'ydir','reverse');
%     legend('RT','MT','Location','southeast');
    

% subplot(4,5,13)
subplot(4,2,[6 8])
    barh(categorical(myalphascales), [rwdpropsyoung(2,:)' rwdpropsyoung(1,:)'] ,'stacked');
    ylabel(sprintf('Reducing Reward\nScale'));
    xlabel(sprintf(['Proportion of savings\n%dJ - %dJ'], rwdhigh, rwdlow));

beautifyfig;

%% Combined figure using old/young data

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
    h(:,i+1)=plot(myalphas,squeeze(old(:,3,alphascaleind,probscaleind,:)).*unit,...
        'Color', 'r', 'LineWidth',1);
    h(2,i+1).LineStyle = '--';
    legend(h(1,[1 3 5 6]),'Young 0.8', 'Young 1.0', 'Young 1.2', 'Old')
    
    patch([myalphas(ptchidx) fliplr(myalphas(ptchidx))], [1.5.*ones(size(myalphas(ptchidx))).*unit -0.5.*ones(size(myalphas(ptchidx))).*unit],...
    ptchcol, 'FaceAlpha',0.3, 'EdgeColor','none','HandleVisibility','off')
    xlabel('Reward (alpha,J)')
    xticks([40 60 80 100])
    yticks([300 600 900 1200 1500])
    ylabel('Duration (ms)')
    set(gca,'ylim',[0.3 1.2].*unit,'xlim',[29 65])
    hold off
    
nexttile(4,[4,3])
    for i = 1:5
        h(:,1)=plot(myalphas,squeeze(young(:,effscaleind,i,probscaleind,:)).*unit,...
                    'LineWidth', 1);
        hold on
        h(2,1).LineStyle = '--';
        h(1,1).Color = i.*[0 0 1/5];
        h(2,1).Color = i.*[0 0 1/5];
    end
    
    patch([myalphas(ptchidx) fliplr(myalphas(ptchidx))], [1.5.*ones(size(myalphas(ptchidx))).*unit -0.5.*ones(size(myalphas(ptchidx))).*unit],...
        ptchcol, 'FaceAlpha',0.3, 'EdgeColor','none')
    legend('RWDScale = 0.8')
    xlabel('Reward (alpha,J)')
    ylabel('Duration (ms)')
    xticks([40 60 80 100])
    yticks([300 600 900 1200 1500])
    set(gca,'ylim',[0.3 1.2].*unit,'xlim',[29 65])
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






