function [J] = utilityfun(x,alpha)
% utilityfun function 
%   Calculates utility for use in comparisons, etc. 
%   Vectorized, unlike mymtrtfun.m which is used in optimization procedure
%   and is solving for MTs and RTs
%   Inputs:
%       - x: nx2 matrix of times. column 1 interpreted as MT, column 2 as RT
%       - alpha: Scalar or nx1 vector
%   Outputs:
%       - J[:,1]: col 1 is the subjective value of reward
%       - J[:,2]: col 2 is the discounted effort
%       - J[:,3]: col 3 is the total utility
%  Dependencies:
%       - global params must be defined beforehand
%       - indcludes a, b, ao, c0, c1, c0t0, c1t0, alphascale,
%       probscale, and effscale

global param

%alpha, c0, c1, c0to, c1to, ao, a, b
Er=@(Tr,a,b) (a.*Tr + b./(Tr)); %reaching cost, basically normal effort cost (no distance, mass, etc.)
Eo=@(To,ao)(ao.*To); %sitting cost, ao is sitting cost, To is sitting
Pr=@(Tr,c0,c1) (1./(1 + exp(-c0 - c1.*Tr))); %Probability of reward as a function of reach time (tr)
Prto =@(To,c0to,c1to) (1./(1 + exp(-c0to - c1to.*To))); %Probability of reward as a function of reacton time (to)

% check dims, make sure its 2 columns
    if size(x,2) == 2
        nreach = size(x,1);
    else 
        x = x';
        nreach = size(x,1);
    end

    if size(alpha,2) ~= 1
        alpha = alpha';
    end
    
    % Initialize storage
    J = zeros([nreach,3]);
    Tr = x(:,1); % MTs
    To = x(:,2); % RTs
    
    % reward rate
    J(:,1) = (param.myalphascale .* alpha .* Pr(Tr, param.myprobscale.*param.myc0, param.myc1) .* Prto(To, param.myc0to, param.myc1to))./(To + Tr);
    
    % effort (W)
    J(:,2) = -1.*(Eo(To, param.myao) + Er(Tr, param.myeffscale.*param.mya, param.myeffscale.*param.myb))./(To + Tr);
    
    % total utility
    J(:,3) = (param.myalphascale .* alpha .* Pr(Tr, param.myprobscale.*param.myc0, param.myc1) .* Prto(To, param.myc0to, param.myc1to)...
        - Eo(To, param.myao) - Er(Tr, param.myeffscale.*param.mya, param.myeffscale.*param.myb))...
        ./(To + Tr);

end

