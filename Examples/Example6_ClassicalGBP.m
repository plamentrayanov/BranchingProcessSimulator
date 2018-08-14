%% EXAMPLE 6: Classical General Branching Process (GBP), also called Crump-Mode-Jagers branching process
% the classical GBP considers time invariant birth and death laws and a population starting from a single individual 
% aged 0 at time 0. This is a toy example that cosinders a human population of only the women. 
% For a more real-life example, you need to extract the life length distribution for men and women, the age specific
% probabilities of giving birth and to use two-type GBP in random environment. Modelling those probabilities of birth 
% and death is done by the classical demography. See Applied Mathematical Demography by Keyfitz and Caswell or the more
% modern methodology on the Human Mortality Database website (https://www.mortality.org/Public/Docs/MethodsProtocol.pdf)

addpath('../')  % adds the function BranchingProcessSimulator to the path
sim_num=1000;   % number of simulations to perform
T=250;          % the horizon is 20 time units and we simulate the branching process in [0, T].
h=0.5;          % the time step is 1 unit of time
omega=min(120, T);        % the maximum age can not be more than the horizon or 120 years
S=(1-normcdf(0:h:omega,76, 10)')./(1-normcdf(0,76, 10));    % the example considers trimmed normal life length
S(end)=0;           % the chance of surviving after age omega is zero

m=0.7;      % the average number of girls a woman has during her life
% less than 1 makes the process subcritical

mu_women_pdf=normpdf(0:h:omega,28, 5)';     % the average age of giving birth is 28 with a standart deviation 5 years
mu_women_pdf([1:12/h, 50/h:end])=0;         % women of age less than 12 ot greater than 50 cannot give birth
mu=m*mu_women_pdf/(sum(mu_women_pdf)*h);    % rescale the trincated distribution so that the total number of children is m

U=1;    % no types - no mutations 

H=[0, 1]';      % women give birth to only 1 child when birth happens
% there is no immigration in the example

Z_0=1;     % the GBP starts from 1 women aged 0 at time 0 

% we can also calculate the age structure of the population which is of interest in this case
% we can use it to calculate the percentage of people on working age, for example
[Z, Z_types, Z_ages] = BranchingProcessSimulator(sim_num, T, h, S, H, U, Z_0, mu);
[Z_mean, Z_lower, Z_upper, Z_median]=confInterval(Z, 0.10);

Z_working=squeeze(sum(Z_ages(:, (18/h):(63/h), 1, :),2));  % population on working age
[Z_working_mean, Z_working_lower, Z_working_upper, Z_working_median]=confInterval(Z_working, 0.10);

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
print('./figures/Example6_fig1', '-dpng', '-r0')

%% number of people on working age
figure('visible','on', 'Units','pixels','OuterPosition',[0 0 1280 1024]);
set(gca,'FontSize',16)
hold on
plot(0:h:T, Z_working', 'Color', [0.7, 0, 0,0.05]);
h_mean=plot(0:h:T, Z_working_mean, 'Color', [0, 0, 0, 0.5], 'LineWidth', line_wd);
h_median=plot(0:h:T, Z_working_median, '--', 'Color', [0, 0, 0, 0.5], 'LineWidth', line_wd);
h_CI=plot(0:h:T, Z_working_lower, '--', 'Color', [0,155/255,1,1], 'LineWidth', line_wd);
plot(0:h:T, Z_working_upper, '--', 'Color', [0,155/255,1,1], 'LineWidth', line_wd);
h_sims=plot(0, 0,'-', 'Color', [0.7, 0, 0], 'LineWidth', line_wd);
legend([h_sims(1), h_mean, h_median, h_CI], 'Simulations', 'Mean', 'Median', '90% conf. interval', 'Location', 'NorthWest')
ylabel('Working Population Count')
xlabel('Time')
print('./figures/Example6_fig2', '-dpng', '-r0')

% save
save(strcat('Example6_', num2str(sim_num)))
