clear all 
% close all

cd("C:\Users\rjc5t\Documents\Neuromechanics\DATA\UtilityModel\rtmtutility\")
load('simresults/simresults_young10.mat'); 
young = mysols;

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

RT = zeros(5);
MT = zeros(5);
RTprops = zeros(5);

for k=1:length(myalphascales)
    for j=1:length(myeffscales)
        % young data
        curvesyoung=squeeze(young(:,j,k,probscaleind,:));
        curvesyoung_mt(:,j)=curvesyoung(:,1);
        curvesyoung_rt(:,j)=curvesyoung(:,2);
        
        savings = ((curvesyoung_rt(alphahigh,j) - curvesyoung_rt(alphalow,j))+(curvesyoung_mt(alphahigh,j) - curvesyoung_mt(alphalow,j))); % total time savings
        %deltasyoung(1,j) = (curvesyoung_mt(alphahigh,j) - curvesyoung_mt(alphalow,j));
        %deltasyoung(2,j) = (curvesyoung_rt(alphahigh,j) - curvesyoung_rt(alphalow,j));
        %propsyoung(1,j) = (curvesyoung_mt(alphahigh,j) - curvesyoung_mt(alphalow,j))/savings; % MT proportion of time savings
        %propsyoung(2,j) = (curvesyoung_rt(alphahigh,j) - curvesyoung_rt(alphalow,j))/savings; % RT proportion of time savings
        RT(j,k) = (curvesyoung_rt(alphahigh,j) - curvesyoung_rt(alphalow,j));
        MT(j,k) = (curvesyoung_mt(alphahigh,j) - curvesyoung_mt(alphalow,j));
        RTprops(j,k) = (curvesyoung_rt(alphahigh,j) - curvesyoung_rt(alphalow,j))/savings;
    end
end

effscalemat = repmat(myeffscales',1,5);
alphascalemat = repmat(myalphascales,5,1);
C = effscalemat.*alphascalemat;

figure;
surf(effscalemat, alphascalemat, RT, C, 'FaceAlpha', 0.5);
xlabel('effort scale'); ylabel('reward scale'); zlabel('\DeltaRT');

figure;
surf(effscalemat, alphascalemat, MT, C, 'FaceAlpha', 0.5);
xlabel('effort scale'); ylabel('reward scale'); zlabel('\DeltaMT');

figure;
surf(effscalemat, alphascalemat, RTprops, C, 'FaceAlpha', 0.5);
zlim([0 1]);
xlabel('effort scale'); ylabel('reward scale'); zlabel('RT proportion');
colorbar('Ticks',[0.7 1.05 1.4],'TickLabels',{'Low J' 'Moderate J' 'High J'});
colormap('spring');
