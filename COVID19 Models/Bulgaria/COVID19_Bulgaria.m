%% Modelling COVID-19 with a General Branching Process (GBP), also called Crump-Mode-Jagers branching process
% The example is built for the specifics of the data and the quarantine in Bulgaria.
% It includes 3 scenarios for the future R0

%% load the branching process simulator
addpath('../../')  % adds the function BranchingProcessSimulator to the path

%% read and plot the data for the chosen country
% save the data from: https://opendata.ecdc.europa.eu/covid19/casedistribution/csv
% then write in the file data.csv
[dates, newcases, totalcases] = readCSVData('data.csv', 'Bulgaria');

%% Build a branching model
detection_time = 8;     % assumes 2 days after the symptoms develop (on average) the person is tested for corona
num_days_passed = dates(end)-dates(1) + 1 + detection_time;   % days passed since the infection began = days since first case + detection time needed

sim_num=1000;    % number of simulations to perform
horizon=60;     % horizon in days, after the last available data
days_before_quarantine=25;      % number of days at the beginning of the process that have a larger R0 value due to less restrictions and social distancing
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
S=(1-normcdf(0:h:omega, 25, 5.1274)')./(1-normcdf(0,25, 5.1274));    % trimmed normal life length with average of 14 days until clearance
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

%% model of the Point process mu
% values are taken from the article https://doi.org/10.1101/2020.02.08.20021162
% from onset to isolation 3.83 (var = 5.99)
% incubation 5.8 (std = 2.6)
% average time of infecting other people is 5.8 + 3.83 = 9.63, std = sqrt(2.6^2 + 5.99) = 3.5707

% R0 - average number of infected people from a single person, conditional on that person staying infected
% distribution of the infection days is taken to be gamma-like and the point process integrates to R0:
% Expectation = k * theta = 9.63, Variance = k * theta^2 = 3.5707^2 => theta = 3.5707^2 / 9.63 = 1.3240, k = 9.63/1.3240 = 7.2734
% Here we assume we have 2 periods of the epidemic in Bulgaria - before the people started distancing themselves and after.
% we have a week of adjustment towards the quarantine where R0 transitions from 3.9 to 0.85
R0 = [ones(1,days_before_quarantine/h)*3.9 linspace(3.9, 0.85, 7/h) ones(1,T/h+1-(7+days_before_quarantine+horizon-detection_time)/h)*0.85 ones(1,(horizon-detection_time)/h)*0.85];
mu_covid_pdf=gampdf(0:h:omega, 7.2734, 1.3240)';
mu_matrix=zeros(size(mu_covid_pdf,1), 1, T/h+1);
mu_matrix(:,1,:)=R0.*repmat(mu_covid_pdf/(sum(mu_covid_pdf)*h),1, T/h+1);
mu=@()(mu_matrix);      % the input format for mu is described in the BranchingProcessSimulator.m

% we can also calculate the age structure of the population which is of interest in this case
% we can use it to calculate the percentage of people on working age, for example
% to get the age structure use: [ActiveCases, ActiveCasesByType, ActiveCasesByAge, TotalCases, TotalCasesByTypes] = BranchingProcessSimulator(sim_num, T, h, S, H, U, Z_0, mu, 'GetAgeStructure', true);
[ActiveCases, ActiveCasesByType, ~, TotalCases, TotalCasesByTypes] = BranchingProcessSimulator(sim_num, T, h, S, H, U, Z_0, mu, 'GetAgeStructure', false);
NewCases = diff(TotalCases');
CuredCases = TotalCases - ActiveCases;

NewCasesPerDay = squeeze(sum(reshape(NewCases, [2, size(NewCases,1)*h, size(NewCases,2)])));
[NewCasesPerDay_mean, NewCasesPerDay_lower, NewCasesPerDay_upper, NewCasesPerDay_median]=confInterval(NewCasesPerDay', 0.10);

%% number of new cases detected per day, model VS reality
% as the detection of the infection is after 8 days, we have 2 Branching Processes here:
% 1) real process happening in the population right now, but not yet detected
% 2) shifted process by 8 days so we can supperimpose with the real data and see if the model parameters are chosen correctly  
line_wd=2.5;
figure('visible','on', 'Units','pixels','OuterPosition',[0 0 1280 1024]);
set(gca,'FontSize',16)
hold on
h_hist=plot(dates, newcases, 'Color', [0, 0, 0, 0.7],'LineWidth', 1.5);
h_median_supperimposed=plot(dates(1):dates(end)+horizon, NewCasesPerDay_median(detection_time:end-1), '-', 'Color', [1, 0, 0, 0.5], 'LineWidth', line_wd);
h_median_real=plot((dates(1):dates(end)+horizon)-detection_time, NewCasesPerDay_median(detection_time:end-1), ':', 'Color', [0, 0.5, 0, 0.5], 'LineWidth', line_wd);
h_CI=plot((dates(1):dates(end)+horizon), NewCasesPerDay_lower(detection_time:end-1), '--', 'Color', [0,155/255,1,1], 'LineWidth', line_wd);
plot((dates(1):dates(end)+horizon), NewCasesPerDay_upper(detection_time:end-1), '--', 'Color', [0,155/255,1,1], 'LineWidth', line_wd);
legend([h_median_real, h_median_supperimposed, h_CI, h_hist], 'Median (Real Process)', 'Prediction of observed new daily cases', '90% conf. interval', 'Daily New Cases', 'Location', 'NorthWest')
xtickangle(90)
dateaxis('X',2)
ylabel('\bf{Daily New Cases}')
xlabel('\bf{Date}')
title('Daily new cases. Real process and timeshifted daily cases prediction')
print('./Figures/Covid_forecast_newcases', '-dpng', '-r0')
%% total number of cases
TotalCasesPerDay= TotalCases(:, 1:(1/h):T/h);
[TotalCases_mean, TotalCases_lower, TotalCases_upper, TotalCases_median]=confInterval(TotalCasesPerDay, 0.10);

figure('visible','on', 'Units','pixels','OuterPosition',[0 0 1280 1024]);
set(gca,'FontSize',16)
hold on
h_hist=plot(dates, cumsum(newcases), 'Color', [0, 0, 0, 0.7],'LineWidth', 1.5);
h_median_supperimposed=plot(dates(1):dates(end)+horizon, TotalCases_median(detection_time:end-1), '-', 'Color', [1, 0, 0, 0.5], 'LineWidth', line_wd);
h_median_real=plot((dates(1):dates(end)+horizon)-detection_time, TotalCases_median(detection_time:end-1), ':', 'Color', [0, 0.5, 0, 0.5], 'LineWidth', line_wd);
h_CI=plot(dates(1):dates(end)+horizon, TotalCases_lower(detection_time:end-1), '--', 'Color', [0,155/255,1,1], 'LineWidth', line_wd);
plot(dates(1):dates(end)+horizon, TotalCases_upper(detection_time:end-1), '--', 'Color', [0,155/255,1,1], 'LineWidth', line_wd);
legend([h_median_real, h_median_supperimposed, h_CI, h_hist], 'Median (Real Process)', 'Prediction of observed total cases', '90% conf. interval', 'Total Cases', 'Location', 'NorthWest')
dateaxis('X',2)
xtickangle(90)
ylabel('\bf{Total Cases}')
xlabel('\bf{Date}')
print('./Figures/Covid_forecast_total', '-dpng', '-r0')

%% total number of active cases
ActiveCasesDaily= ActiveCases(:, 1:(1/h):T/h);
[ActiveCasesDaily_mean, ActiveCasesDaily_lower, ActiveCasesDaily_upper, ActiveCasesDaily_median]=confInterval(ActiveCasesDaily, 0.10);

figure('visible','on', 'Units','pixels','OuterPosition',[0 0 1280 1024]);
set(gca,'FontSize',16)
hold on
h_median_supperimposed=plot(dates(1):dates(end)+horizon, ActiveCasesDaily_median(detection_time:end-1), '-', 'Color', [1, 0, 0, 0.5], 'LineWidth', line_wd);
h_CI=plot(dates(1):dates(end)+horizon, ActiveCasesDaily_lower(detection_time:end-1), '--', 'Color', [0,155/255,1,1], 'LineWidth', line_wd);
plot(dates(1):dates(end)+horizon, ActiveCasesDaily_upper(detection_time:end-1), '--', 'Color', [0,155/255,1,1], 'LineWidth', line_wd);
legend([h_median_supperimposed, h_CI], 'Prediction of observed active cases', '90% conf. interval', 'Location', 'NorthWest')
dateaxis('X',2)
xtickangle(90)
ylabel('\bf{Active Cases}')
xlabel('\bf{Date}')
print('./Figures/Covid_forecast_active', '-dpng', '-r0')

%% OPTIMISTIC SCENARIO model of the Point process mu, Optimistic Scenario, R0=0.8, from now on
% values are taken from the article https://doi.org/10.1101/2020.02.08.20021162
% from onset to isolation 3.83 (var = 5.99)
% incubation 5.8 (std = 2.6)
% average time of infecting other people is 5.8 + 3.83 = 9.63, std = sqrt(2.6^2 + 5.99) = 3.5707

% R0 - average number of infected people from a single person, conditional on that person staying infected
% distribution of the infection days is taken to be gamma-like and the point process integrates to R0:
% Expectation = k * theta = 9.63, Variance = k * theta^2 = 3.5707^2 => theta = 3.5707^2 / 9.63 = 1.3240, k = 9.63/1.3240 = 7.2734
% Here we assume we have 2 periods of the epidemic in Bulgaria - before the people started distancing themselves and after.
% we have a week of adjustment towards the quarantine where R0 transitions from 3.9 to 0.85
R0 = [ones(1,days_before_quarantine/h)*3.9 linspace(3.9, 0.85, 7/h) ones(1,T/h+1-(7+days_before_quarantine+horizon-detection_time)/h)*0.85 ones(1,(horizon-detection_time)/h)*0.7];
mu_covid_pdf=gampdf(0:h:omega, 7.2734, 1.3240)';
mu_matrix=zeros(size(mu_covid_pdf,1), 1, T/h+1);
mu_matrix(:,1,:)=R0.*repmat(mu_covid_pdf/(sum(mu_covid_pdf)*h),1, T/h+1);
mu=@()(mu_matrix);      % the input format for mu is described in the BranchingProcessSimulator.m

% we can also calculate the age structure of the population which is of interest in this case
% we can use it to calculate the percentage of people on working age, for example
[ActiveCases, ActiveCasesByType, ~, TotalCases, TotalCasesByTypes] = BranchingProcessSimulator(sim_num, T, h, S, H, U, Z_0, mu, 'GetAgeStructure', false);
NewCases = diff(TotalCases');
CuredCases = TotalCases - ActiveCases;

NewCasesPerDay = squeeze(sum(reshape(NewCases, [2, size(NewCases,1)*h, size(NewCases,2)])));
[NewCasesPerDay_mean, NewCasesPerDay_lower, NewCasesPerDay_upper, NewCasesPerDay_median]=confInterval(NewCasesPerDay', 0.10);

%% OPTIMISTIC SCENARIO number of new cases detected per day, model VS reality
% as the detection of the infection is after 8 days, we have 2 Branching Processes here:
% 1) real process happening in the population right now, but not yet detected
% 2) shifted process by 8 days so we can supperimpose with the real data and see if the model parameters are chosen correctly  
line_wd=2.5;
figure('visible','on', 'Units','pixels','OuterPosition',[0 0 1280 1024]);
set(gca,'FontSize',16)
hold on
h_hist=plot(dates, newcases, 'Color', [0, 0, 0, 0.7],'LineWidth', 1.5);
h_median_supperimposed=plot(dates(1):dates(end)+horizon, NewCasesPerDay_median(detection_time:end-1), '-', 'Color', [1, 0, 0, 0.5], 'LineWidth', line_wd);
h_median_real=plot((dates(1):dates(end)+horizon)-detection_time, NewCasesPerDay_median(detection_time:end-1), ':', 'Color', [0, 0.5, 0, 0.5], 'LineWidth', line_wd);
h_CI=plot((dates(1):dates(end)+horizon), NewCasesPerDay_lower(detection_time:end-1), '--', 'Color', [0,155/255,1,1], 'LineWidth', line_wd);
plot((dates(1):dates(end)+horizon), NewCasesPerDay_upper(detection_time:end-1), '--', 'Color', [0,155/255,1,1], 'LineWidth', line_wd);
legend([h_median_real, h_median_supperimposed, h_CI, h_hist], 'Median (Real Process)', 'Prediction of observed new daily cases', '90% conf. interval', 'Daily New Cases', 'Location', 'NorthWest')
xtickangle(90)
dateaxis('X',2)
ylabel('\bf{Daily New Cases}')
xlabel('\bf{Date}')
title('Daily new cases. Real process and timeshifted daily cases prediction')
print('./Figures/Covid_forecast_newcases_optimistic', '-dpng', '-r0')
%% OPTIMISTIC SCENARIO total number of cases
TotalCasesPerDay= TotalCases(:, 1:(1/h):T/h);
[TotalCases_mean, TotalCases_lower, TotalCases_upper, TotalCases_median]=confInterval(TotalCasesPerDay, 0.10);

figure('visible','on', 'Units','pixels','OuterPosition',[0 0 1280 1024]);
set(gca,'FontSize',16)
hold on
h_hist=plot(dates, cumsum(newcases), 'Color', [0, 0, 0, 0.7],'LineWidth', 1.5);
h_median_supperimposed=plot(dates(1):dates(end)+horizon, TotalCases_median(detection_time:end-1), '-', 'Color', [1, 0, 0, 0.5], 'LineWidth', line_wd);
h_median_real=plot((dates(1):dates(end)+horizon)-detection_time, TotalCases_median(detection_time:end-1), ':', 'Color', [0, 0.5, 0, 0.5], 'LineWidth', line_wd);
h_CI=plot(dates(1):dates(end)+horizon, TotalCases_lower(detection_time:end-1), '--', 'Color', [0,155/255,1,1], 'LineWidth', line_wd);
plot(dates(1):dates(end)+horizon, TotalCases_upper(detection_time:end-1), '--', 'Color', [0,155/255,1,1], 'LineWidth', line_wd);
legend([h_median_real, h_median_supperimposed, h_CI, h_hist], 'Median (Real Process)', 'Prediction of observed total cases', '90% conf. interval', 'Total Cases', 'Location', 'NorthWest')
dateaxis('X',2)
xtickangle(90)
ylabel('\bf{Total Cases}')
xlabel('\bf{Date}')
print('./Figures/Covid_forecast_total_optimistic', '-dpng', '-r0')

%% OPTIMISTIC SCENARIO total number of active cases
ActiveCasesDaily= ActiveCases(:, 1:(1/h):T/h);
[ActiveCasesDaily_mean, ActiveCasesDaily_lower, ActiveCasesDaily_upper, ActiveCasesDaily_median]=confInterval(ActiveCasesDaily, 0.10);

figure('visible','on', 'Units','pixels','OuterPosition',[0 0 1280 1024]);
set(gca,'FontSize',16)
hold on
h_median_supperimposed=plot(dates(1):dates(end)+horizon, ActiveCasesDaily_median(detection_time:end-1), '-', 'Color', [1, 0, 0, 0.5], 'LineWidth', line_wd);
h_CI=plot(dates(1):dates(end)+horizon, ActiveCasesDaily_lower(detection_time:end-1), '--', 'Color', [0,155/255,1,1], 'LineWidth', line_wd);
plot(dates(1):dates(end)+horizon, ActiveCasesDaily_upper(detection_time:end-1), '--', 'Color', [0,155/255,1,1], 'LineWidth', line_wd);
legend([h_median_supperimposed, h_CI], 'Prediction of observed active cases', '90% conf. interval', 'Location', 'NorthWest')
dateaxis('X',2)
xtickangle(90)
ylabel('\bf{Active Cases}')
xlabel('\bf{Date}')
print('./Figures/Covid_forecast_active_optimistic', '-dpng', '-r0')

%% PESSIMISTIC SCENARIO model of the Point process mu, Pessimistic Scenario, R0=1.5, from now on
% values are taken from the article https://doi.org/10.1101/2020.02.08.20021162
% from onset to isolation 3.83 (var = 5.99)
% incubation 5.8 (std = 2.6)
% average time of infecting other people is 5.8 + 3.83 = 9.63, std = sqrt(2.6^2 + 5.99) = 3.5707

% R0 - average number of infected people from a single person, conditional on that person staying infected
% distribution of the infection days is taken to be gamma-like and the point process integrates to R0:
% Expectation = k * theta = 9.63, Variance = k * theta^2 = 3.5707^2 => theta = 3.5707^2 / 9.63 = 1.3240, k = 9.63/1.3240 = 7.2734
% Here we assume we have 2 periods of the epidemic in Bulgaria - before the people started distancing themselves and after.
% we have a week of adjustment towards the quarantine where R0 transitions from 3.9 to 0.85
R0 = [ones(1,days_before_quarantine/h)*3.9 linspace(3.9, 0.85, 7/h) ones(1,T/h+1-(7+days_before_quarantine+horizon-detection_time)/h)*0.85 ones(1,(horizon-detection_time)/h)*1.3];
mu_covid_pdf=gampdf(0:h:omega, 7.2734, 1.3240)';
mu_matrix=zeros(size(mu_covid_pdf,1), 1, T/h+1);
mu_matrix(:,1,:)=R0.*repmat(mu_covid_pdf/(sum(mu_covid_pdf)*h),1, T/h+1);
mu=@()(mu_matrix);      % the input format for mu is described in the BranchingProcessSimulator.m

% we can also calculate the age structure of the population which is of interest in this case
% we can use it to calculate the percentage of people on working age, for example
[ActiveCases, ActiveCasesByType, ~, TotalCases, TotalCasesByTypes] = BranchingProcessSimulator(sim_num, T, h, S, H, U, Z_0, mu, 'GetAgeStructure', false);
NewCases = diff(TotalCases');
CuredCases = TotalCases - ActiveCases;

NewCasesPerDay = squeeze(sum(reshape(NewCases, [2, size(NewCases,1)*h, size(NewCases,2)])));
[NewCasesPerDay_mean, NewCasesPerDay_lower, NewCasesPerDay_upper, NewCasesPerDay_median]=confInterval(NewCasesPerDay', 0.10);

%% PESSIMISTIC SCENARIO number of new cases detected per day, model VS reality
% as the detection of the infection is after 8 days, we have 2 Branching Processes here:
% 1) real process happening in the population right now, but not yet detected
% 2) shifted process by 8 days so we can supperimpose with the real data and see if the model parameters are chosen correctly  
line_wd=2.5;
figure('visible','on', 'Units','pixels','OuterPosition',[0 0 1280 1024]);
set(gca,'FontSize',16)
hold on
h_hist=plot(dates, newcases, 'Color', [0, 0, 0, 0.7],'LineWidth', 1.5);
h_median_supperimposed=plot(dates(1):dates(end)+horizon, NewCasesPerDay_median(detection_time:end-1), '-', 'Color', [1, 0, 0, 0.5], 'LineWidth', line_wd);
h_median_real=plot((dates(1):dates(end)+horizon)-detection_time, NewCasesPerDay_median(detection_time:end-1), ':', 'Color', [0, 0.5, 0, 0.5], 'LineWidth', line_wd);
h_CI=plot((dates(1):dates(end)+horizon), NewCasesPerDay_lower(detection_time:end-1), '--', 'Color', [0,155/255,1,1], 'LineWidth', line_wd);
plot((dates(1):dates(end)+horizon), NewCasesPerDay_upper(detection_time:end-1), '--', 'Color', [0,155/255,1,1], 'LineWidth', line_wd);
legend([h_median_real, h_median_supperimposed, h_CI, h_hist], 'Median (Real Process)', 'Prediction of observed new daily cases', '90% conf. interval', 'Daily New Cases', 'Location', 'NorthWest')
xtickangle(90)
dateaxis('X',2)
ylabel('\bf{Daily New Cases}')
xlabel('\bf{Date}')
title('Daily new cases. Real process and timeshifted daily cases prediction')
print('./Figures/Covid_forecast_newcases_pessimistic', '-dpng', '-r0')
%% PESSIMISTIC SCENARIO total number of cases
TotalCasesPerDay= TotalCases(:, 1:(1/h):T/h);
[TotalCases_mean, TotalCases_lower, TotalCases_upper, TotalCases_median]=confInterval(TotalCasesPerDay, 0.10);

figure('visible','on', 'Units','pixels','OuterPosition',[0 0 1280 1024]);
set(gca,'FontSize',16)
hold on
h_hist=plot(dates, cumsum(newcases), 'Color', [0, 0, 0, 0.7],'LineWidth', 1.5);
h_median_supperimposed=plot(dates(1):dates(end)+horizon, TotalCases_median(detection_time:end-1), '-', 'Color', [1, 0, 0, 0.5], 'LineWidth', line_wd);
h_median_real=plot((dates(1):dates(end)+horizon)-detection_time, TotalCases_median(detection_time:end-1), ':', 'Color', [0, 0.5, 0, 0.5], 'LineWidth', line_wd);
h_CI=plot(dates(1):dates(end)+horizon, TotalCases_lower(detection_time:end-1), '--', 'Color', [0,155/255,1,1], 'LineWidth', line_wd);
plot(dates(1):dates(end)+horizon, TotalCases_upper(detection_time:end-1), '--', 'Color', [0,155/255,1,1], 'LineWidth', line_wd);
legend([h_median_real, h_median_supperimposed, h_CI, h_hist], 'Median (Real Process)', 'Prediction of observed total cases', '90% conf. interval', 'Total Cases', 'Location', 'NorthWest')
dateaxis('X',2)
xtickangle(90)
ylabel('\bf{Total Cases}')
xlabel('\bf{Date}')
print('./Figures/Covid_forecast_total_pessimistic', '-dpng', '-r0')

%% OPTIMISTIC SCENARIO total number of active cases
ActiveCasesDaily= ActiveCases(:, 1:(1/h):T/h);
[ActiveCasesDaily_mean, ActiveCasesDaily_lower, ActiveCasesDaily_upper, ActiveCasesDaily_median]=confInterval(ActiveCasesDaily, 0.10);

figure('visible','on', 'Units','pixels','OuterPosition',[0 0 1280 1024]);
set(gca,'FontSize',16)
hold on
h_median_supperimposed=plot(dates(1):dates(end)+horizon, ActiveCasesDaily_median(detection_time:end-1), '-', 'Color', [1, 0, 0, 0.5], 'LineWidth', line_wd);
h_CI=plot(dates(1):dates(end)+horizon, ActiveCasesDaily_lower(detection_time:end-1), '--', 'Color', [0,155/255,1,1], 'LineWidth', line_wd);
plot(dates(1):dates(end)+horizon, ActiveCasesDaily_upper(detection_time:end-1), '--', 'Color', [0,155/255,1,1], 'LineWidth', line_wd);
legend([h_median_supperimposed, h_CI], 'Prediction of observed active cases', '90% conf. interval', 'Location', 'NorthWest')
dateaxis('X',2)
xtickangle(90)
ylabel('\bf{Active Cases}')
xlabel('\bf{Date}')
print('./Figures/Covid_forecast_active_pessimistic', '-dpng', '-r0')

%% WORST CASE SCENARIO model of the Point process mu, Pessimistic Scenario, R0=1.5, from now on
% values are taken from the article https://doi.org/10.1101/2020.02.08.20021162
% from onset to isolation 3.83 (var = 5.99)
% incubation 5.8 (std = 2.6)
% average time of infecting other people is 5.8 + 3.83 = 9.63, std = sqrt(2.6^2 + 5.99) = 3.5707

% R0 - average number of infected people from a single person, conditional on that person staying infected
% distribution of the infection days is taken to be gamma-like and the point process integrates to R0:
% Expectation = k * theta = 9.63, Variance = k * theta^2 = 3.5707^2 => theta = 3.5707^2 / 9.63 = 1.3240, k = 9.63/1.3240 = 7.2734
% Here we assume we have 2 periods of the epidemic in Bulgaria - before the people started distancing themselves and after.
% we have a week of adjustment towards the quarantine where R0 transitions from 3.9 to 0.85
R0 = [ones(1,days_before_quarantine/h)*3.9 linspace(3.9, 0.85, 7/h) ones(1,T/h+1-(7+days_before_quarantine+horizon-detection_time)/h)*0.85 ones(1,(horizon-detection_time)/h)*3];
mu_covid_pdf=gampdf(0:h:omega, 7.2734, 1.3240)';
mu_matrix=zeros(size(mu_covid_pdf,1), 1, T/h+1);
mu_matrix(:,1,:)=R0.*repmat(mu_covid_pdf/(sum(mu_covid_pdf)*h),1, T/h+1);
mu=@()(mu_matrix);      % the input format for mu is described in the BranchingProcessSimulator.m

% we can also calculate the age structure of the population which is of interest in this case
% we can use it to calculate the percentage of people on working age, for example
[ActiveCases, ActiveCasesByType, ~, TotalCases, TotalCasesByTypes] = BranchingProcessSimulator(sim_num, T, h, S, H, U, Z_0, mu, 'GetAgeStructure', false);
NewCases = diff(TotalCases');
CuredCases = TotalCases - ActiveCases;

NewCasesPerDay = squeeze(sum(reshape(NewCases, [2, size(NewCases,1)*h, size(NewCases,2)])));
[NewCasesPerDay_mean, NewCasesPerDay_lower, NewCasesPerDay_upper, NewCasesPerDay_median]=confInterval(NewCasesPerDay', 0.10);

%% WORST CASE SCENARIO number of new cases detected per day, model VS reality
% as the detection of the infection is after 8 days, we have 2 Branching Processes here:
% 1) real process happening in the population right now, but not yet detected
% 2) shifted process by 8 days so we can supperimpose with the real data and see if the model parameters are chosen correctly  
line_wd=2.5;
figure('visible','on', 'Units','pixels','OuterPosition',[0 0 1280 1024]);
set(gca,'FontSize',16)
hold on
h_hist=plot(dates, newcases, 'Color', [0, 0, 0, 0.7],'LineWidth', 1.5);
h_median_supperimposed=plot(dates(1):dates(end)+horizon, NewCasesPerDay_median(detection_time:end-1), '-', 'Color', [1, 0, 0, 0.5], 'LineWidth', line_wd);
h_median_real=plot((dates(1):dates(end)+horizon)-detection_time, NewCasesPerDay_median(detection_time:end-1), ':', 'Color', [0, 0.5, 0, 0.5], 'LineWidth', line_wd);
h_CI=plot((dates(1):dates(end)+horizon), NewCasesPerDay_lower(detection_time:end-1), '--', 'Color', [0,155/255,1,1], 'LineWidth', line_wd);
plot((dates(1):dates(end)+horizon), NewCasesPerDay_upper(detection_time:end-1), '--', 'Color', [0,155/255,1,1], 'LineWidth', line_wd);
legend([h_median_real, h_median_supperimposed, h_CI, h_hist], 'Median (Real Process)', 'Prediction of observed new daily cases', '90% conf. interval', 'Daily New Cases', 'Location', 'NorthWest')
xtickangle(90)
dateaxis('X',2)
ylabel('\bf{Daily New Cases}')
xlabel('\bf{Date}')
title('Daily new cases. Real process and timeshifted daily cases prediction')
print('./Figures/Covid_forecast_newcases_worstcase', '-dpng', '-r0')
%% WORST CASE SCENARIO total number of cases
TotalCasesPerDay= TotalCases(:, 1:(1/h):T/h);
[TotalCases_mean, TotalCases_lower, TotalCases_upper, TotalCases_median]=confInterval(TotalCasesPerDay, 0.10);

figure('visible','on', 'Units','pixels','OuterPosition',[0 0 1280 1024]);
set(gca,'FontSize',16)
hold on
h_hist=plot(dates, cumsum(newcases), 'Color', [0, 0, 0, 0.7],'LineWidth', 1.5);
h_median_supperimposed=plot(dates(1):dates(end)+horizon, TotalCases_median(detection_time:end-1), '-', 'Color', [1, 0, 0, 0.5], 'LineWidth', line_wd);
h_median_real=plot((dates(1):dates(end)+horizon)-detection_time, TotalCases_median(detection_time:end-1), ':', 'Color', [0, 0.5, 0, 0.5], 'LineWidth', line_wd);
h_CI=plot(dates(1):dates(end)+horizon, TotalCases_lower(detection_time:end-1), '--', 'Color', [0,155/255,1,1], 'LineWidth', line_wd);
plot(dates(1):dates(end)+horizon, TotalCases_upper(detection_time:end-1), '--', 'Color', [0,155/255,1,1], 'LineWidth', line_wd);
legend([h_median_real, h_median_supperimposed, h_CI, h_hist], 'Median (Real Process)', 'Prediction of observed total cases', '90% conf. interval', 'Total Cases', 'Location', 'NorthWest')
dateaxis('X',2)
xtickangle(90)
ylabel('\bf{Total Cases}')
xlabel('\bf{Date}')
print('./Figures/Covid_forecast_total_worstcase', '-dpng', '-r0')
%% WORST CASE SCENARIO total number of active cases
ActiveCasesDaily= ActiveCases(:, 1:(1/h):T/h);
[ActiveCasesDaily_mean, ActiveCasesDaily_lower, ActiveCasesDaily_upper, ActiveCasesDaily_median]=confInterval(ActiveCasesDaily, 0.10);

figure('visible','on', 'Units','pixels','OuterPosition',[0 0 1280 1024]);
set(gca,'FontSize',16)
hold on
h_median_supperimposed=plot(dates(1):dates(end)+horizon, ActiveCasesDaily_median(detection_time:end-1), '-', 'Color', [1, 0, 0, 0.5], 'LineWidth', line_wd);
h_CI=plot(dates(1):dates(end)+horizon, ActiveCasesDaily_lower(detection_time:end-1), '--', 'Color', [0,155/255,1,1], 'LineWidth', line_wd);
plot(dates(1):dates(end)+horizon, ActiveCasesDaily_upper(detection_time:end-1), '--', 'Color', [0,155/255,1,1], 'LineWidth', line_wd);
legend([h_median_supperimposed, h_CI], 'Prediction of observed active cases', '90% conf. interval', 'Location', 'NorthWest')
dateaxis('X',2)
xtickangle(90)
ylabel('\bf{Active Cases}')
xlabel('\bf{Date}')
print('./Figures/Covid_forecast_active_worstcase', '-dpng', '-r0')

%% Shows the simulations with confidence intervals
% to plot the Active cases use:

%%% MU, POINT PROCESS
line_wd=3.5;
figure('visible','on', 'Units','pixels','OuterPosition',[0 0 1280 1024]);
set(gca,'FontSize',16)
hold on
plot(0:h:omega, squeeze(mu()), 'LineWidth', line_wd);
ylabel('Point process densities')
xlabel('Days')
title('Probability for spreading the infection, by days')
print('./Figures/Covid_PointProcessDensity', '-dpng', '-r0')
%%%

%%% Infected period distribution
line_wd=3.5;
figure('visible','on', 'Units','pixels','OuterPosition',[0 0 1280 1024]);
set(gca,'FontSize',16)
hold on
plot(0:h:omega, [diff(1-S); 0]./h, 'LineWidth', line_wd);
ylabel('Infected period PDF')
xlabel('Days')
title('Probability for getting cured of the infection, by days')
print('./Figures/Covid_InfectedPeriodPDF', '-dpng', '-r0')
%%%


%% save the results
save(strcat('COVID19_', num2str(sim_num)))

