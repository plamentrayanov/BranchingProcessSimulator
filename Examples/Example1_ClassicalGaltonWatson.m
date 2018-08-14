%% EXAMPLE 1: The Classical Galton-Watson process
addpath('../')  % adds the function BranchingProcessSimulator to the path
sim_num=100000; % number of simulations to perform
T=50;   % the horizon is 50 time units and we simulate the branching process in [0, T].
h=1;    % the time step is 1 unit of time
S=[1,0]';  % survivability function defines the particles to die exactly after 1 unit of time
Z_0=1;     % the initial population starts with 1 particle

U=1;   % we have only 1 type, so the mutation probabilities matrix is just 1
H=[0.55,0,0.45]';   % using the classical notation, the offspring p.g.f. is defined as f(s)=0.55 + 0*s^1 + 0.45*s^2
m=(0:2)*H;          % the mean number of offspring
% for classical GW process, there is no immigration and no point process

% we are not interested in the age structure, as it is trivial in this case
[Z, Z_types] = BranchingProcessSimulator(sim_num, T, h, S, H, U, Z_0);
[Z_mean_100, Z_lower_100, Z_upper_100, Z_median_100]=confInterval(Z(1:100,:), 0.10);
[Z_mean_1000, Z_lower_1000, Z_upper_1000, Z_median_1000]=confInterval(Z(1:1000,:), 0.10);
[Z_mean_10000, Z_lower_10000, Z_upper_10000, Z_median_10000]=confInterval(Z(1:10000,:), 0.10);
[Z_mean_100000, Z_lower_100000, Z_upper_100000, Z_median_100000]=confInterval(Z(1:100000,:), 0.10);

% it is known from the classical theory that the mean population of such process at time t is m^t
Z_mean_true=m.^(0:h:T);

%% Shows the convergence of the simulated mean toward the true mean
figure('visible','on', 'Units','pixels','OuterPosition',[0 0 1280 1024]);
set(gca,'FontSize',16)
line_wd=2.5;
hold on
%plot(0:h:T, Z(1:1000,:)', 'Color', [0.7, 0, 0,0.05]);
h_mean_100=plot(0:h:T, Z_mean_100, 'Color', [1, 0, 0, 0.5], 'LineWidth', line_wd);
h_mean_1000=plot(0:h:T, Z_mean_1000, 'Color', [0, 0.7, 0, 0.5], 'LineWidth', line_wd);
h_mean_10000=plot(0:h:T, Z_mean_10000, 'Color', [0, 0, 1, 0.5], 'LineWidth', line_wd);
h_mean_100000=plot(0:h:T, Z_mean_100000, 'Color', [247/250, 150/250, 70/250, 0.5], 'LineWidth', line_wd);
h_mean_true=plot(0:h:T, Z_mean_true, ':', 'Color', [0, 0, 0, 0.5], 'LineWidth', line_wd);
legend([h_mean_100, h_mean_1000, h_mean_10000, h_mean_100000, h_mean_true], 'Sim Mean (100)', 'Sim Mean (1000)', ...
    'Sim Mean (10000)', 'Sim Mean (100000)', 'True mean (m^t)', 'Location', 'NorthEast')
ylabel('Expected Population Count')
xlabel('Time')
print('./figures/Example1_fig1', '-dpng', '-r0')

%% draw confidence intervals and median
figure('visible','on', 'Units','pixels','OuterPosition',[0 0 1280 1024], 'PaperPositionMode','auto');
set(gca,'FontSize',16)
hold on
plot(0:h:T, Z(1:1000,:)', 'Color', [0.7, 0, 0,0.05]);
h_mean=plot(0:h:T, Z_mean_100000, 'Color', [0, 0, 0, 0.5], 'LineWidth', line_wd);
h_median=plot(0:h:T, Z_median_100000, 'x', 'Color', [0, 0, 0, 0.5], 'LineWidth', line_wd);
h_CI=plot(0:h:T, Z_lower_100000, '--', 'Color', [0,155/255,1,1], 'LineWidth', line_wd);
plot(0:h:T, Z_upper_100000, '--', 'Color', [0,155/255,1,1], 'LineWidth', line_wd);
h_sims=plot(0, Z(1,1),'-', 'Color', [0.7, 0, 0], 'LineWidth', line_wd);
legend([h_sims(1), h_mean, h_median, h_CI], 'Simulations', 'Mean', 'Median', '90% conf. interval', 'Location', 'NorthWest')
ylabel('Total Population Count')
xlabel('Time')
print('./figures/Example1_fig2', '-dpng', '-r0')

% save
save(strcat('Example1_', num2str(sim_num)))
