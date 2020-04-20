%% Modelling COVID-19 with a General Branching Process (GBP), also called Crump-Mode-Jagers branching process
% The example is built for the specifics of the data and the quarantine in Bulgaria.
% It includes 3 scenarios for the future R0

%% load the branching process simulator
addpath('../../')       % adds the function BranchingProcessSimulator to the path
addpath('../')          % adds the functions confInterval and readCSVData to the path

%% read and plot the data for the chosen country
% save the data from: https://opendata.ecdc.europa.eu/covid19/casedistribution/csv
% then write in the file data.csv
[dates, newcases_hist, totalcases_hist] = readCSVData('../data.csv', 'Sweden');

% for Sweden, the first cases seem controlled and should be excluded as a beginning of the epidemics
% the epidemics starts on :
ep_start = find(newcases_hist==0, 1, 'last')+1;
dates=dates(ep_start:end);
newcases_hist=newcases_hist(ep_start:end);
totalcases_hist=totalcases_hist(ep_start:end);
%% Build a branching model
detection_time = 8;     % assumes 2 days after the symptoms develop (on average) the person is tested for corona
num_days_passed = dates(end)-dates(1) + 1 + detection_time;   % days passed since the infection began = days since first case + detection time needed

sim_num=1000;    % number of simulations to perform
horizon=90;     % horizon in days, after the last available data
%days_before_quarantine=20;      % number of days at the beginning of the process that have a larger R0 value due to less restrictions and social distancing
T=num_days_passed+horizon;          % the simulation period for the process, from the first infected to horizon 
h=0.5;          % the time step is 0.5 day
omega=60;           % no one is infected for more than omega days, technical parameter that does not need to be accurate, smaller -> faster speed, 
                    % but it needs to be at least as large that the P(being infected after omega) is very close to zero.

%% model of S, survivability function, NOT to be mistaken with the percentage of people who survive the virus!
% S is virus specific and represents the survivability function of the virus inside the body of a person
% Here it is assumed normally distributed where the values mu and sigma are taken from the article:
% "Feasability of controlling 2019-nCoV outbreaks by isolation of cases and contacts", J. Hellewell, S. Abbot, et al. 
% https://doi.org/10.1101/2020.02.08.20021162
%
% average time for recovery in the article is 5.8 + 9.1 = 14.9, std = sqrt(2.6^2 + 19.53) = 5.1274
% however the data in bulgaria (recovery rates and active cases) seems to show that people need more time to clear the virus out of their body
% even after they are already healthy and recovered. Here a value of 25 days is used for the average days to clear the virus.
% As worldometers.info writes, the active cases and recovories depend highly on the way we count a recovery.
% In Bulgaria it is counted after a person tests negative twice for the virus.
S=(1-normcdf(0:h:omega, 30, 5.1274)')./(1-normcdf(0,30, 5.1274));    % trimmed normal life length with average of 14 days until clearance
S(end)=0;           % the chance of surviving after age omega is zero
% If you want to try exponential distribution instead, use S=(1-expcdf(0:h:omega, 14.9)')./normcdf(0,omega, 14.9);   % exponential with average of 14 days

%% other parameters for the simulator
U=1;    % no types - no mutations, there is a possibility to include a second, quarantined type here with a probability to become quarantined

H=[0, 1]';      % people infect only 1 other person when infection happens, no multiple infections

% Im=[];
% there is no immigration in the example, but it could be included as a random variable representing the number of individuals returning to Bulgaria
% which are infected with the virus

% the initial number of infections is uniformly distributed with respect to "days passed since infection" (i.e. "age" of the infection)
% those infections have uniformly distributed "ages" on every call of the Z_0 function
age_struct_pdf=unifpdf(0:h:omega, 0, 5.8)';
age_struct_prob=age_struct_pdf./(sum(age_struct_pdf));

Z_0=@()(mnrnd(5, age_struct_prob)');

%% MAIN SCENARIO, model of the Point process mu
% values are taken from the article https://doi.org/10.1101/2020.02.08.20021162
% from onset to isolation 3.83 (var = 5.99)
% incubation 5.8 (std = 2.6)
% average time of infecting other people is 5.8 + 3.83 = 9.63, std = sqrt(2.6^2 + 5.99) = 3.5707

% R0 - average number of infected people from a single person, conditional on that person staying infected
% distribution of the infection days is taken to be gamma-like and the point process integrates to R0:
% Expectation = k * theta = 9.63, Variance = k * theta^2 = 3.5707^2 => theta = 3.5707^2 / 9.63 = 1.3240, k = 9.63/1.3240 = 7.2734
% Here we assume we have 2 periods of the epidemic - before the people started distancing themselves and after.
% R0 is chosen to fit what happened already
R0 = [linspace(7, 5.5, 19/h) linspace(5.5, 2.3, 13/h) linspace(2.3, 1.9, T/h+1-(24+7+25+horizon-detection_time)/h) linspace(1.9, 0.9, 19/h) linspace(0.9, 0.9, 5/h) ones(1,(horizon-detection_time)/h)*0.9];
mu_covid_pdf=gampdf(0:h:omega, 7.2734, 1.3240)';
mu_matrix=zeros(size(mu_covid_pdf,1), 1, T/h+1);
mu_matrix(:,1,:)=R0.*repmat(mu_covid_pdf/(sum(mu_covid_pdf)*h),1, T/h+1);
mu=@()(mu_matrix);      % the input format for mu is described in the BranchingProcessSimulator.m

% we can also calculate the age structure of the population which is of interest in this case
% we can use it to calculate the percentage of people on working age, for example
% to get the age structure use: [ActiveCases, ActiveCasesByType, ActiveCasesByAge, TotalCases, TotalCasesByTypes] = BranchingProcessSimulator(sim_num, T, h, S, H, U, Z_0, mu, 'GetAgeStructure', true);
[ActiveCases, ~, ~, TotalCases, ~] = BranchingProcessSimulator(sim_num, T, h, S, H, U, Z_0, mu, 'GetAgeStructure', false);
NewCases = diff(TotalCases');
% CuredCases = TotalCases - ActiveCases;

% converts from h period of time to daily period of time
NewCasesDaily = squeeze(sum(reshape(NewCases, [2, size(NewCases,1)*h, size(NewCases,2)])))';
TotalCasesDaily= TotalCases(:, 1:(1/h):T/h);    
ActiveCasesDaily= ActiveCases(:, 1:(1/h):T/h);

% build and save plots
buildPlots(NewCasesDaily, TotalCasesDaily, ActiveCasesDaily, newcases_hist, dates, horizon, detection_time, ...
            'MainScenario', 'No change in R_0 (no change in measures)');

%% OPTIMISTIC SCENARIO model of the Point process mu, Optimistic Scenario, R0 changes from 0.8 to 0.6, from now on
R0 = [linspace(7, 5.5, 19/h) linspace(5.5, 2.3, 13/h) linspace(2.3, 1.9, T/h+1-(24+7+25+horizon-detection_time)/h) linspace(1.9, 0.9, 19/h) linspace(0.9, 0.9, 5/h) ones(1,(horizon-detection_time)/h)*0.7];
mu_covid_pdf=gampdf(0:h:omega, 7.2734, 1.3240)';
mu_matrix=zeros(size(mu_covid_pdf,1), 1, T/h+1);
mu_matrix(:,1,:)=R0.*repmat(mu_covid_pdf/(sum(mu_covid_pdf)*h),1, T/h+1);
mu_optimistic=@()(mu_matrix);      % the input format for mu is described in the BranchingProcessSimulator.m

[ActiveCases_optimistic, ~, ~, TotalCases_optimistic, ~] = BranchingProcessSimulator(sim_num, T, h, S, H, U, Z_0, mu_optimistic, 'GetAgeStructure', false);
NewCases_optimistic = diff(TotalCases_optimistic');
% CuredCases = TotalCases - ActiveCases;

% converts from h period of time to daily period of time
NewCasesDaily_optimistic = squeeze(sum(reshape(NewCases_optimistic, [2, size(NewCases_optimistic,1)*h, size(NewCases_optimistic,2)])))';
TotalCasesDaily_optimistic= TotalCases(:, 1:(1/h):T/h);    
ActiveCasesDaily_optimistic= ActiveCases(:, 1:(1/h):T/h);

% build and save plots
buildPlots(NewCasesDaily_optimistic, TotalCasesDaily_optimistic, ActiveCasesDaily_optimistic, newcases_hist, dates, horizon, detection_time, ...
            'OptimisticScenario', 'Decline in R_0 (better results from measures)');

%% PESSIMISTIC SCENARIO model of the Point process mu, Pessimistic Scenario, R0 changes from 0.8 to 1.2, from now on
R0 = [linspace(7, 5.5, 19/h) linspace(5.5, 2.3, 13/h) linspace(2.3, 1.9, T/h+1-(24+7+25+horizon-detection_time)/h) linspace(1.9, 0.9, 19/h) linspace(0.9, 0.9, 5/h) ones(1,(horizon-detection_time)/h)*1.3];
mu_covid_pdf=gampdf(0:h:omega, 7.2734, 1.3240)';
mu_matrix=zeros(size(mu_covid_pdf,1), 1, T/h+1);
mu_matrix(:,1,:)=R0.*repmat(mu_covid_pdf/(sum(mu_covid_pdf)*h),1, T/h+1);
mu_pessimistic=@()(mu_matrix);      % the input format for mu is described in the BranchingProcessSimulator.m

[ActiveCases_pessimistic, ~, ~, TotalCases_pessimistic, ~] = BranchingProcessSimulator(sim_num, T, h, S, H, U, Z_0, mu_pessimistic, 'GetAgeStructure', false);
NewCases_pessimistic = diff(TotalCases_pessimistic');
% CuredCases = TotalCases - ActiveCases;

% converts from h period of time to daily period of time
NewCasesDaily_pessimistic = squeeze(sum(reshape(NewCases_pessimistic, [2, size(NewCases_pessimistic,1)*h, size(NewCases_pessimistic,2)])))';
TotalCasesDaily_pessimistic= TotalCases(:, 1:(1/h):T/h);    
ActiveCasesDaily_pessimistic= ActiveCases(:, 1:(1/h):T/h);

% build and save plots
buildPlots(NewCasesDaily_pessimistic, TotalCasesDaily_pessimistic, ActiveCasesDaily_pessimistic, newcases_hist, dates, horizon, detection_time, ...
            'PessimisticScenario', 'Increase in R_0 (worse results from measures)');

%% Comparison of scenarios
[NewCasesDaily_mean, NewCasesDaily_lower, NewCasesDaily_upper, NewCasesDaily_median]=confInterval(NewCasesDaily, 0.10);
[NewCasesDaily_optimisic_mean, NewCasesDaily_optimisic_lower, NewCasesDaily_optimisic_upper, NewCasesDaily_optimisic_median]=confInterval(NewCasesDaily_optimistic, 0.10);
[NewCasesDaily_pessimistic_mean, NewCasesDaily_pessimistic_lower, NewCasesDaily_pessimistic_upper, NewCasesDaily_pessimistic_median]=confInterval(NewCasesDaily_pessimistic, 0.10);

line_wd=2.5;
figure('visible','on', 'Units','pixels','OuterPosition',[0 0 1280 1024]);
set(gca,'FontSize',16)
hold on
h_hist=plot(dates, newcases_hist, 'Color', [0, 0, 0, 0.7],'LineWidth', 1.5);
h_main=plot(dates(1):dates(end)+horizon, NewCasesDaily_median(detection_time:end-1), '-', 'Color', [0, 0, 0.5, 0.5], 'LineWidth', line_wd);
h_optimistic=plot(dates(end):dates(end)+horizon, NewCasesDaily_optimisic_median(end-horizon-1:end-1), '-', 'Color', [0, 0.5, 0, 0.5], 'LineWidth', line_wd);
h_pessimistic=plot(dates(end):dates(end)+horizon, NewCasesDaily_pessimistic_median(end-horizon-1:end-1), '-', 'Color', [0.5, 0, 0, 0.5], 'LineWidth', line_wd);
h_CI=plot(dates(1):dates(end)+horizon, NewCasesDaily_lower(detection_time:end-1), '--', 'Color', [0,155/255,1,1], 'LineWidth', line_wd);
plot(dates(1):dates(end)+horizon, NewCasesDaily_upper(detection_time:end-1), '--', 'Color', [0,155/255,1,1], 'LineWidth', line_wd);
legend([h_main, h_optimistic, h_pessimistic, h_hist, h_CI], 'Main Scenario', 'Optimistic Scenario', 'Pessimistic Scenario', 'Observed new daily cases', '90% conf. interval', 'Location', 'NorthWest')
xtickangle(90)
x_ticks = xticks;
xticks(x_ticks(1):7:x_ticks(end))
dateaxis('X',2)
ylabel('\bf{Total Cases (Observed)}')
xlabel('\bf{Date}')
title('Forecasts of observed New Daily Cases (by scenario)')
print(strcat('./Figures/forecast_newcases_by_scenario'), '-dpng', '-r0')

%% Shows the input to the branching process
%%% MU, POINT PROCESS
line_wd=2.5;
figure('visible','on', 'Units','pixels','OuterPosition',[0 0 1280 1024]);
set(gca,'FontSize',16)
hold on
mu_forecast = squeeze(mu());
plot(0:h:omega, mu_forecast(:, end), 'Color', [0.5,0,0,0.5], 'LineWidth', line_wd);
ylabel('Point process densities')
xlabel('Days')
title('Probability for spreading the infection, by days')
print('./Figures/PointProcessDensity', '-dpng', '-r0')

%%% Infected period distribution
line_wd=2.5;
figure('visible','on', 'Units','pixels','OuterPosition',[0 0 1280 1024]);
set(gca,'FontSize',16)
hold on
plot(0:h:omega, [diff(1-S); 0]./h, 'Color', [0,0,0.5,0.5], 'LineWidth', line_wd);
ylabel('Infected period PDF')
xlabel('Days')
title('Probability for getting cleared of the virus, by days')
print('./Figures/InfectedPeriodPDF', '-dpng', '-r0')

%% save the results
% save(strcat('Italy_', num2str(sim_num)))

