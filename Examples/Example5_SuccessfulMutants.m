%% EXAMPLE 5: Multi-type Bellman-Harris with successful mutants
% This example simulates the model considered in article:
% "Branching processes in continuous time as models of mutations: Computational approaches and algorithms
% by Slavtchova-Bojkova, M., Trayanov, P., Dimitrov, S. (2017). 
% Computational Statistics and Data Analysis, http://dx.doi.org/10.1016/j.csda.2016.12.013.
%
% The article model description in short:
% type 0 particles - supercritical, do not mutate to other types once they occur
% type 1 particles - subcritical, with probability of mutation to type 0
% The population starts from 1 particle of type 1.
% "Successful mutant" - a particle of type 0 that starts a type 0 branching, that never extincts
% "Unsuccessful mutant" - a particle of type 0 that starts a type 0 branching, that extincts
% The object of interest in the article is to calculate the distribution T of the first arrival of a "successful mutant".
%
% This script shows 2 different ways to simulate the model and how the distribution of T can be found by simulations.
% 

addpath('../')  % adds the function BranchingProcessSimulator to the path
sim_num=1000;   % number of simulations to perform
T=150;          % the horizon, we simulate the branching process in [0, T].
h=0.1;          % the time step is 0.1 
omega=T;        % maximum possible age in [0, T]
S=[[1, 0, zeros(1, omega/h-1)]', (1-normcdf(0:h:omega,10, 2.5)')./(1-normcdf(0,10, 2.5))];  % it is explained in the 
% article that the survivability function of type 0 cells does not matter to the distribution of T, so any distribution 
% will be ok. Type 1 has truncated normal distribution of the lifelength.

Z_0=[0 1]';     % we start from a single type 1 particle

q_0=0.3;        % type 0 is defined to have probability of extinction 0.3 
H=[q_0/(1+q_0), 0, 1/(1+q_0); 1-0.375, 0, 0.375]';  % type 1 is subcritical, doomed to extinction
% the p.g.f. of type 0 particles is defined in such a way that the probability for extinction is exactly q_0!
% Again, any p.g.f. for which this is true, will do!

u=0.2;      % mutation probability from type 1 to type 0
U=[1 0; u, 1-u]';    % type 0 is not mutating

[~, Z_types] = BranchingProcessSimulator(sim_num, T, h, S, H, U, Z_0);

% calculating however the distribution P(T>t) from this simulation is troublesome, as we do not get as output 
% which mutants are successful and which are not. We can approximate that if a simulated type 0 population is non-zero 
% at the end of the horizon, then the population is never going to extinct and the mutant that started it is successful

% get the finite times of occurence of successful mutants
time_of_occurance=cellfun(@(x) find(x, 1, 'first'), num2cell(squeeze(Z_types(Z_types(:,1,end)>0,1,:)), 2))*h;
% then add the cases in which T=inf:
time_of_occurance=[time_of_occurance; inf(sim_num-length(time_of_occurance), 1)];

% calculate the empirical CDF
[F_occurance1, X_occurance1]=ecdf(time_of_occurance);

%% An alternative simulation scheme with better approximation
% instead of considering two-type branching process, we can define an analogous scheme using 3 types
% type 2: the subcritical process (the previous type 1)
% type 1: unsuccessful mutants
% type 0: successful mutants

S=[[1, 0, zeros(1, omega/h-1)]', [1, 0, zeros(1, omega/h-1)]', ...  % type 0 and 1 die shortly after they are born
    (1-normcdf(0:h:omega,10, 2.5)')./(1-normcdf(0,10, 2.5))];       % type 2 has truncated normal life length

Z_0=[0 0 1]';   % again, we start from the subcritical process
q_0=0.3;        % the same extinction probability

% type 0 and 1 have no offspring, type 2 has the same offspring distribution as before
H=[1,0,0; 1,0,0; 1-0.375, 0, 0.375]';
u=0.2;  % mutation probability

% instead of using the previous mutation probability matrix, we redefine it: type 0 and type 1 do not mutate 
% (actually they also produce no offspring so it does not matter what you define for them here)
% type 2 produces a mutant with probability u (type 0 or type 1) and this mutant is successful if it does not lead to 
% extinction (which happens with probability q_0). So the mutation probability from type 2 to type 0 is (1-q_0)*u 
% and the mutation probability to type 1 is q_0*u.
U=[1 0 0; 0 1 0; (1-q_0)*u, q_0*u, 1-u]';

% lets simulate now the 3-type branching process
[~, Z_types] = BranchingProcessSimulator(sim_num, T, h, S, H, U, Z_0);

% get the finite times of occurence of successful mutants
time_of_occurance=cellfun(@(x) find(x, 1, 'first'), num2cell(squeeze(Z_types(sum(Z_types(:, 1, :), 3)>0,1,:)), 2))*h;
% then add the cases in which T=inf:
time_of_occurance=[time_of_occurance; inf(sim_num-length(time_of_occurance), 1)];

% calculate the empirical CDF
[F_occurance2, X_occurance2]=ecdf(time_of_occurance);

%% plot P(T>t) calculated from both simulation schemes
line_wd=2.5;
figure('visible','on', 'Units','pixels','OuterPosition',[0 0 1280 1024]);
hold on
set(gca,'FontSize',16)
plot(X_occurance1, 1-F_occurance1, 'Color', [0, 0, 0, 0.5], 'LineWidth', line_wd)
plot(X_occurance2, 1-F_occurance2, 'Color', [1, 0, 0, 0.5], 'LineWidth', line_wd)
ylabel('P(T>t)')
xlabel('Time')
print('./figures/Example5_fig1', '-dpng', '-r0')

% The second scheme is actually closer to the true answer, as P(T=inf)=0.815.
% The 2 schemes however should produce similar results as sim_num gets larger

%% save
save(strcat('Example5_', num2str(sim_num)))
