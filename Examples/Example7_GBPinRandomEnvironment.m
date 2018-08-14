%% EXAMPLE 7: Multi-type General Branching Process (GBP) in random environment
% the example demonstrates the input format for constant, varying and random environment:
%
% U is taken to be constant non-random
%
% H is taken to be varying through time, but non-random, i.e. the same for every draw
%
% Z_0 is taken to be random, drawn when the branching begins, ages are uniformly distributed
%
% S is taken to be random, stationary stochastic process
%
% mu is taken to be random, non-stationary stochastic process
%
% Im is taken to be random, homogeneous across time, with rate of 0.10 persons per year, unformly aged
%
% Note that all of the simulation paths for U, H, mu, S, Im and Z_0 are independent of eachother and of the branching 
% process itself.

addpath('../')  % adds the function BranchingProcessSimulator to the path
sim_num=1000;   % number of simulations to perform
T=50;          % the horizon is 20 time units and we simulate the branching process in [0, T].
h=0.5;          % the time step is 1 unit of time
omega=120;        % the maximum age

U=[1 0; 0, 1]';    % no types - no mutations 

% define H to be different 2D matrix for each t (varying, non-random environment):
p=linspace(0, 0.2, T/h);
H=zeros(3,2,T/h);
for t=1:size(H, 3)  % go through time, from 0 to T
    H(:, :, t)=[1, 0, 0; 0, 1-p(t), p(t)]';
end

% the initial number of individuals is at least 1, binomially distributed with parameters n=10, p=0.5
% those individuals have uniformly distributed ages on every call of the Z_0 function
Z_0=@()([zeros(omega/h+1, 1), mnrnd(1+binornd(100,0.5), ones(1,omega/h+1)./(omega/h+1))']);

% we can also calculate the age structure of the population which is of interest in this case
% we can use it to calculate the percentage of people on working age, for example
[Z, Z_types, Z_ages] = BranchingProcessSimulator(sim_num, T, h, @()(draw_S(T, h, omega)), H, U, Z_0, @()(draw_mu(T, h, omega)), @()(draw_Im(T, h, omega)));
[Z_mean, Z_lower, Z_upper, Z_median]=confInterval(Z, 0.10);

% men get pension after age 65, women after age 63
Z_working=squeeze(sum(Z_ages(:, (18/h):(63/h), 2, :),2) + sum(Z_ages(:, (18/h):(65/h), 1, :),2));  % population on working age
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
print('./figures/Example7_fig1', '-dpng', '-r0')

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
print('./figures/Example7_fig2', '-dpng', '-r0')

% save
save(strcat('Example7_', num2str(sim_num)))

%% functions definitions 
% model S to has static noise around it
function S=draw_S(T, h, omega)
% generates truncated normal survival functions for men and women with means and STDs that are uniformly distributed
% on every draw. The survival functions are modelled to be subjected to random noise. S is modelled as a stationary 
% stochastic process.
% OMEGA is the highest age in the life table.

S=zeros(omega/h+1, 2, T/h+1);
    for t=1:T/h
        mean_w=unifrnd(70, 74);     % mean life length for women is uniformly distributed between 70 and 74
        std_w=unifrnd(8,12);        % the STD is also subject to noise
        mean_m=unifrnd(74, 78);     % mean life length for men is uniformly distributed between 74 and 78
        std_m=unifrnd(9,11);
        % generate a truncated normal distribution for the life lengths using the simulated means and STDs
        S(:,:,t)=[(1-normcdf(0:h:omega, mean_w, std_w)')./(1-normcdf(0,mean_w, std_w)), (1-normcdf(0:h:omega,mean_m, std_m)')./(1-normcdf(0,mean_m, std_m))];
        S(end,:,t)=0;   % the survival probability beyond age omega is 0
    end
end

% model S to has a trend and noise in the mean, i.e. non-stationary stochastic process
function mu=draw_mu(T, h, omega)
% generate the whole simulation path for mu, on every call
mu_mean=linspace(28, 35, T/h);      % the average age of giving birth is moving from 28 to 35 years-old in [0, T]
mu=zeros(omega/h+1, 2, T/h);        % 
for t=1:T/h
    m=unifrnd(1.5, 2);      % the average number of children is subject to uniform noise, but between 1.5 and 2
    mu_women_pdf=normpdf(0:h:omega, mu_mean(t), 5)';    % the density of mu is taken normal
    mu_women_pdf([1:12/h, 50/h:end])=0;     % women do not give birth before age 12 and after age 50
    mu(:,:,t)=[zeros(size(mu_women_pdf)), m*mu_women_pdf/(sum(mu_women_pdf)*h)];    % men do not give birth and
    % the pdf for women is rescaled, so the women have an average number of m children through their lives
end
end

% model the immigration
function Im=draw_Im(T, h, omega)
% generate the whole simulation path of immigrations on every call
Im=zeros(omega/h+1, 2, T/h+1);    % no immigration at time 0
for t=1:T/h+1
    % no immigration for men, uniform distributed ages for the women immigrants
    % the immigration rate is 1 person per 10 years
    Im(:,:,t)=[zeros(omega/h+1, 1), mnrnd(binornd(1, h*0.1), ones(1,omega/h+1)./(omega/h+1))'];
end
end



