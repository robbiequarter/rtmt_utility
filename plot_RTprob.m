% plot RT probability functions

% global parameters needed
    global param
    param.myc0to = -1; % Used in RT function
    param.myc1to = 10; % Used in RT function
    param.myc0to2 = -6; % Shallower - used in new RT function
    param.myc1to2 = 20; % Shallower - used in new RT function
%     param.myc0to2 = -8.75; % Steeper - used in new RT function
%     param.myc1to2 = 35; % Steeper - used in new RT function
    param.ntarg = ceil((2*pi*0.1)/0.016);

% Scaled up RT probability function
    Pr_oldyoung = @(Tr) (1./(1 + exp(-(param.myc0to) - (param.myc1to).*Tr)));

    % Alternative RT function
    k00 = 20; 
    x00 = 0.30;
    Pr_test2 = @(Tr) ((1-(1/param.ntarg))./(1 + exp(-(param.myc0to2) - (param.myc1to2).*Tr))) + (1/param.ntarg); % n targets
    Pr_test3 = @(Tr) (0.75./(1 + exp(-(param.myc0to2) - (param.myc1to2).*Tr))) + 0.25; % 4 targets

% Plot it
    figure
    fplot(Pr_oldyoung,'Color','k', 'LineWidth', 1, 'DisplayName', 'Orig. RT function')
    hold on
    fplot(Pr_test2, 'Color', [.7 .7 .7], 'LineWidth', 1, 'DisplayName', 'N targ. RT function')
    fplot(Pr_test3, 'Color', [.8 .8 .8], 'LineWidth', 1, 'DisplayName', '4 targ. RT function')
    xline(0,'k--','HandleVisibility','off'); 
    yline(.25,'k--','HandleVisibility','off');
    legend('show', 'Location','southeast');
    xlim([-0.5 1.0]);
    xlabel('Reaction time (s)'); ylabel('Probability of success');
    hold off