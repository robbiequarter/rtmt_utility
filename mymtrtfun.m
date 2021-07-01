function J = mymtrtfun(x)
% mymtrtfun function 
%   Calculates utility with a given To (rxn time) and Tr (mvt duration)
%   Used in conjunction with mtrt models loop 
%   Being used in fminsearch()

global param

%alpha, c0, c1, c0to, c1to, ao, a, b
Er=@(Tr,a,b) (a*Tr + b./(Tr)); %reaching cost, basically normal effort cost (no distance, mass, etc.)
Eo=@(To,ao)(ao*To); %sitting cost, ao is sitting cost, To is sitting
% Pr=@(Tr,c0,c1) 1; %Probability of reward as a function of reach time (tr)
Pr=@(Tr,c0,c1) (1./(1 + exp(-c0 - c1*Tr))); %Probability of reward as a function of reach time (tr)
Prto =@(To,c0to,c1to) (1./(1 + exp(-c0to - c1to*To))); %(*Probability of reward as a function of reacton time (to)

% New RT function:
% Prto =@(To,c0to,c1to)(0.75./(1 + exp(-(c0to) - (c1to).*To)))+0.25; %(*Probability of reward as a function of reacton time (to)

Tr=x(1); % MT
To=x(2); % RT

% %make negative since we are really after maximum, but will be minimizing this function
% Effort scale on a and b
J=-(param.myalphascale * param.myalpha * Pr(Tr,param.myprobscale*param.myc0,param.myc1)*Prto(To,param.myc0to,param.myc1to)...
    - Eo(To, param.myao) - Er(Tr, param.myeffscale * param.mya, param.myeffscale * param.myb))/(To + Tr);
% Effort scale only on b
% J=-(param.myalphascale * param.myalpha * Pr(Tr,param.myprobscale*param.myc0,param.myc1)*Prto(To,param.myc0to,param.myc1to)...
%     - Eo(To, param.myao) - Er(Tr, param.mya, param.myeffscale * param.myb))/(To + Tr);

end

