%% EXAMPLE 2: The Classical Galton-Watson process
addpath('../')  % adds the function BranchingProcessSimulator to the path
sim_num=1000;   % number of simulations to perform
T=20;   % the horizon is 20 time units and we simulate the branching process in [0, T].
h=1;    % the time step is 1 unit of time
S=[ones(1,3); zeros(1,3)];  % all 3 types of the particles die exactly after 1 unit of time
Z_0=[0,500,0]';     % the initial population starts with 500 particles of type 2

U=[0.9,0.05,0.05; 0.05,0.95,0; 0.05,0,0.95]';   % mutation probabilities matrix
H=[1/3,0,2/3; 3/4,0,1/4; 1/2,0,1/2]';   % type 1 is supercritical, type 2 - subcritical, type 3 - critical
% there is no immigration and no point process in this example

% we are not interested in the age structure, as it is trivial in this case
[Z, Z_types] = BranchingProcessSimulator(sim_num, T, h, S, H, U, Z_0);
[Z_mean, Z_lower, Z_upper, Z_median]=confInterval(Z, 0.10);

%% Shows the simulations with confidence intervals
line_wd=2.5;
figure('visible','on', 'Units','pixels','OuterPosition',[0 0 1280 1024]);
set(gca,'FontSize',16)
hold on
plot(0:h:T, Z', 'Color', [0.7, 0, 0,0.05]);
h_mean=plot(0:h:T, Z_mean, 'Color', [0, 0, 0, 0.5], 'LineWidth', line_wd);
h_median=plot(0:h:T, Z_median, '--', 'Color', [0, 0, 0, 0.5], 'LineWidth', line_wd);
h_CI=plot(0:h:T, Z_lower, '--', 'Color', [0,155/255,1,1], 'LineWidth', line_wd);
plot(0:h:T, Z_upper, '--', 'Color', [0,155/255,1,1], 'LineWidth', line_wd);
h_sims=plot(0, Z(1,1),'-', 'Color', [0.7, 0, 0], 'LineWidth', line_wd);
legend([h_sims(1), h_mean, h_median, h_CI], 'Simulations', 'Mean', 'Median', '90% conf. interval', 'Location', 'NorthWest')
ylabel('Total Population Count')
xlabel('Time')
print('./figures/Example2_fig1', '-dpng', '-r0')

%% draw the average number of particles by type
figure('visible','on', 'Units','pixels','OuterPosition',[0 0 1280 1024]);
set(gca,'FontSize',16)
hold on
area(0:h:T, squeeze(mean(Z_types,1))')
colormap([0:(1/(size(Z_types,2)-1)):1; zeros(2, size(Z_types,2))]')
legend(strcat({'Type '}, cellstr(num2str((1:size(Z_types,2))'))), 'Location', 'NorthWest')
ylabel('Population Count by Type')
xlabel('Time')
print('./figures/Example2_fig2', '-dpng', '-r0')

% save
save(strcat('Example2_', num2str(sim_num)))
