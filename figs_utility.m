% plotting and comparing changes in effort/reward states and overall
% utility

clear all 
close all

cd("C:\Users\rjc5t\Documents\Neuromechanics\DATA\UtilityModel\rtmtutility\")
load('simresults/simresults_young13.mat'); 
young = mysols;

%% Params
global param 
param.myc0 = -5; % accuracy parameters; shifts logistic to the right with scaling
param.myc1 = 10; % accuracy parameters; lower value reduces steepness of logistic curve
param.mya = 77; % effort offset 
param.myb = 12; %b in metabolic equation (new = 11)
param.myao = 77; % resting rate for younger adults
param.myi=1.23;  %exponent on distance
param.myc0to = -1; % Used in optimization RT accuracy
param.myc1to = 10; % Used in optimization RT accuracy

%range of alpha values
myalphas = [20:200]';

%range of SCALING on effort, alpha, and probability terms
param.myalphascale = 1;
param.myeffscale = 1;
param.myprobscale = 1;
myeffscales = 0.8:0.1:1.2;
myalphascales= 0.8:0.1:1.2;
myprobscales = 0.8:0.1:1.2;
alphascaleind=find(myalphascales==1);
effscaleind=find(myeffscales==1);
probscaleind= find(myprobscales==1);

% Alphas of certain indices
nrwdlow = 35;
rwdlow = 40;
nrwdlowidx = find(myalphas==nrwdlow);
rwdlowidx = find(myalphas==rwdlow);

nrwdhigh = 195;
rwdhigh = 200;
nrwdhighidx = find(myalphas==nrwdhigh);
rwdhighidx = find(myalphas==rwdhigh);

alphas = [nrwdlow; rwdlow; nrwdhigh; rwdhigh];

%% Loop through and calculate J's at certain alphas
effJ = [];
rwdJ = [];
probJ = [];

for j = 1:5
    effcurvesyoung = squeeze(young(:,j,alphascaleind,probscaleind,:));
    effcurvesyoung_mt(:,j) = effcurvesyoung(:,1);
    effcurvesyoung_rt(:,j) = effcurvesyoung(:,2);
    Jtemp = utilityfun(effcurvesyoung([nrwdlowidx, rwdlowidx, nrwdhighidx, rwdhighidx], :), alphas);
    effJ(:,j) = Jtemp(:,3);
    effE(:,j) = Jtemp(:,2);
    effR(:,j) = Jtemp(:,1);
    
    rwdcurvesyoung = squeeze(young(:,effscaleind,j,probscaleind,:));
    rwdcurvesyoung_mt(:,j) = rwdcurvesyoung(:,1);
    rwdcurvesyoung_rt(:,j) = rwdcurvesyoung(:,2);
    Jtemp = utilityfun(rwdcurvesyoung([nrwdlowidx, rwdlowidx, nrwdhighidx, rwdhighidx], :), alphas);
    rwdJ(:,j) = Jtemp(:,3);
    rwdE(:,j) = Jtemp(:,2);
    rwdR(:,j) = Jtemp(:,1);
    
    probcurvesyoung = squeeze(young(:,effscaleind,alphascaleind,j,:));
    probcurvesyoung_mt(:,j) = probcurvesyoung(:,1);
    probcurvesyoung_rt(:,j) = probcurvesyoung(:,2);    
    Jtemp = utilityfun(probcurvesyoung([nrwdlowidx, rwdlowidx, nrwdhighidx, rwdhighidx], :), alphas);
    probJ(:,j) = Jtemp(:,3);
    probE(:,j) = Jtemp(:,2);
    probR(:,j) = Jtemp(:,1);
end

%% Plot J
figure
subplot(1,3,1)
hold on
    plot(myeffscales, effR, 'r.-'); 
    plot(myeffscales, rwdR, 'g.-');
    plot(myeffscales, probR, 'b.-');
    xlabel('Scaling Factor'); ylabel('Subjective Value'); 

subplot(1,3,2)
hold on
    p(:,1) = plot(myeffscales, effE, 'r.-');
    p(:,2) = plot(myeffscales, rwdE, 'g.-');
    p(:,3) = plot(myeffscales, probE, 'b.-');
    xlabel('Scaling Factor'); ylabel('Effort');
    legend([p(1,1), p(1,2), p(1,3)],{'Effort', 'RWD', 'Prob'},'Location','southeast')
    
subplot(1,3,3)
hold on
    plot(myeffscales, effJ, 'r.-') 
    plot(myeffscales, rwdJ, 'g.-')
    plot(myeffscales, probJ, 'b.-')
    xlabel('Scaling Factor'); ylabel('Utility');

%% Bar plot of delta Js

scale = 5;
deltaJ = [effJ(2,scale) - effJ(1,scale),... % lower RWD - lower NRWD
    rwdJ(2,1) - rwdJ(1,1),...
    probJ(2,scale) - probJ(1,scale);
    effJ(4,scale) - effJ(3,scale),... % higher RWD - higher NRWD
    rwdJ(4,1) - rwdJ(3,1),...
    probJ(4,scale) - probJ(3,scale)];

figure
bar(categorical({sprintf('%dJ - %dJ',rwdlow, nrwdlow),sprintf('%dJ - %dJ',rwdhigh, nrwdhigh)}),deltaJ);
legend('Effort Scale = 1.2', 'RWD Scale = 0.8','Prob Scale = 1.2','Location','northoutside');
xlabel('\DeltaRWD'); ylabel('\DeltaJ');
