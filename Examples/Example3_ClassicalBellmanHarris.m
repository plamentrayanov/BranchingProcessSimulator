%% EXAMPLE 3: The Classical Galton-Watson process
addpath('../')      % adds the function BranchingProcessSimulator to the path
sim_num=10000;      % number of simulations to perform
T=20;       % the horizon is 20 time units and we simulate the branching process in [0, T].
h=0.1;      % the time step is 0.1
omega=T;    % oldest age in the population
mean_lifelength=3;      % defines the average lifelength
S=1-expcdf(0:h:omega, mean_lifelength)';    % the lifelength is exponentially distributed
Z_0=1;    % the initial population starts with 1 particle of age 0 at time 0

U=1;    % no mutations and types
H=[3/4,0,1/4]';   % the process is subcritical
% there is no immigration and no point process in this example
m=(0:2)*H;      % mean number of children
Z_mean_true=exp((m-1)*(0:h:T)./mean_lifelength);     % the theoretical solution

% we are not interested in the age structure, as it is trivial in this case
[Z, Z_types] = BranchingProcessSimulator(sim_num, T, h, S, H, U, Z_0);
[Z_mean_100, Z_lower_100, Z_upper_100, Z_median_100]=confInterval(Z(1:100,:), 0.10);
[Z_mean_1000, Z_lower_1000, Z_upper_1000, Z_median_1000]=confInterval(Z(1:1000,:), 0.10);
[Z_mean_10000, Z_lower_10000, Z_upper_10000, Z_median_10000]=confInterval(Z(1:10000,:), 0.10);

%% Shows the convergence of the simulated mean toward the true mean
figure('visible','on', 'Units','pixels','OuterPosition',[0 0 1280 1024]);
set(gca,'FontSize',16)
line_wd=2.5;
hold on
%plot(0:h:T, Z(1:1000,:)', 'Color', [0.7, 0, 0,0.05]);
h_mean_100=plot(0:h:T, Z_mean_100, 'Color', [1, 0, 0, 0.5], 'LineWidth', line_wd);
h_mean_1000=plot(0:h:T, Z_mean_1000, 'Color', [0, 0.7, 0, 0.5], 'LineWidth', line_wd);
h_mean_10000=plot(0:h:T, Z_mean_10000, 'Color', [0, 0, 1, 0.5], 'LineWidth', line_wd);
h_mean_true=plot(0:h:T, Z_mean_true, ':', 'Color', [0, 0, 0, 0.5], 'LineWidth', line_wd);
legend([h_mean_100, h_mean_1000, h_mean_10000, h_mean_true], 'Sim Mean (100)', 'Sim Mean (1000)', ...
    'Sim Mean (10000)', 'True mean', 'Location', 'NorthEast')
ylabel('Expected Population Count')
xlabel('Time')
print('./figures/Example3_fig1', '-dpng', '-r0')

%% draw confidence intervals and median
figure('visible','on', 'Units','pixels','OuterPosition',[0 0 1280 1024], 'PaperPositionMode','auto');
set(gca,'FontSize',16)
hold on
plot(0:h:T, Z(1:1000,:)', 'Color', [0.7, 0, 0,0.05]);
h_mean=plot(0:h:T, Z_mean_10000, 'Color', [0, 0, 0, 0.5], 'LineWidth', line_wd);
h_median=plot(0:h:T, Z_median_10000, ':', 'Color', [0, 0, 0, 0.5], 'LineWidth', line_wd);
h_CI=plot(0:h:T, Z_lower_10000, '--', 'Color', [0,155/255,1,1], 'LineWidth', line_wd);
plot(0:h:T, Z_upper_10000, '--', 'Color', [0,155/255,1,1], 'LineWidth', line_wd);
h_sims=plot(0, Z(1,1),'-', 'Color', [0.7, 0, 0], 'LineWidth', line_wd);
legend([h_sims(1), h_mean, h_median, h_CI], 'Simulations', 'Mean', 'Median', '90% conf. interval', 'Location', 'NorthWest')
ylabel('Total Population Count')
xlabel('Time')
print('./figures/Example3_fig2', '-dpng', '-r0')

% save
save(strcat('Example3_', num2str(sim_num)))
