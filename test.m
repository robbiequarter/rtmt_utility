% young
load('simresults/simresults_young.mat'); 
young = mysols;
% old
load('simresults/simresults_old.mat'); 
old = mysols;

global param
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

%alpha, c0, c1, c0to, c1to, ao, a, b
Er=@(Tr,a,b) (a*Tr + b./(Tr)); %reaching cost, basically normal effort cost (no distance, mass, etc.)
Eo=@(To,ao)(ao*To); %sitting cost, ao is sitting cost, To is sitting
Pr=@(Tr,c0,c1) 1; %Probability of reward as a function of reach time (tr)
Prto =@(To,c0to,c1to) (1./(1 + exp(-c0to - c1to*To))); %(*Probability of reward as a function of reacton time (to)

J = @(alpha,To,Tr,effscale) -(alpha *1*Prto(To,param.myc0to,param.myc1to)...
    - Eo(To, param.myaoyoung) - Er(Tr, effscale*param.myayoung, effscale*param.mybyoung))/(To + Tr);

% Eff Scale = 1
MT1 = squeeze(young(:,3,alphascaleind,probscaleind,1)); % MT
RT1 = squeeze(young(:,3,alphascaleind,probscaleind,2)); % RT
for i = 11:length(myalphas)
    utils(i-10,1) = J(myalphas(i), RT1(i), MT1(i), 1);
end

% Eff Scale = 1.2
MT2 = squeeze(young(:,5,alphascaleind,probscaleind,1)); % MT
RT2 = squeeze(young(:,5,alphascaleind,probscaleind,2)); % RT
for i = 11:length(myalphas)
    utils(i-10,2) = J(myalphas(i), RT2(i), MT2(i), 1.2);
end

%% Plots
subplot(2,2,1); plot(RT1(11:end), utils(:,1), RT2(11:end), utils(:,2));
xlabel('RT (s)'); ylabel('Utility (J/s)'); legend('1.0','1.2'); 
subplot(2,2,3); plot(MT1(11:end), utils(:,1), MT2(11:end), utils(:,2));
xlabel('MT (s)'); ylabel('Utility (J/s)'); legend('1.0','1.2');
subplot(2,2,[2,4]); plot(RT1(11:end)+MT1(11:end), utils(:,1), RT2(11:end)+MT2(11:end), utils(:,2));
xlabel('Total Time (s)'); ylabel('Utility (J/s)'); legend('1.0','1.2'); 

figure
subplot(2,2,1); plot(utils(:,1), RT1(11:end),utils(:,2), RT2(11:end));
xlabel('RT (s)'); ylabel('Utility (J/s)'); legend('1.0','1.2'); 
subplot(2,2,3); plot(utils(:,1), MT1(11:end), utils(:,2), MT2(11:end));
xlabel('MT (s)'); ylabel('Utility (J/s)'); legend('1.0','1.2');
subplot(2,2,[2,4]); plot(utils(:,1), RT1(11:end)+MT1(11:end),utils(:,2), RT2(11:end)+MT2(11:end));
xlabel('Total Time (s)'); ylabel('Utility (J/s)'); legend('1.0','1.2'); 

figure
subplot(2,2,1); plot(myalphas(11:end),RT1(11:end), myalphas(11:end),RT2(11:end));
ylabel('RT (s)'); xlabel('RWD (J)'); legend('1.0','1.2'); 
subplot(2,2,3); plot(myalphas(11:end),MT1(11:end), myalphas(11:end),MT2(11:end));
ylabel('MT (s)'); xlabel('RWD (J)'); legend('1.0','1.2');
subplot(2,2,[2,4]); plot(myalphas(11:end),RT1(11:end)+MT1(11:end), myalphas(11:end), RT2(11:end)+MT2(11:end));
ylabel('Total Time (s)'); xlabel('RWD (J)'); legend('1.0','1.2'); 
