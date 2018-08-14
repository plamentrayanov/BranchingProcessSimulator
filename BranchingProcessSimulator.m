function [Z, Z_types, Z_ages]=BranchingProcessSimulator(sim_num, T, h, draw_S, draw_H, draw_U, draw_Z_0, varargin)
% Written in Matlab
% AUTHOR: Plamen Trayanov
% Faculty of Mathematics and Informatics
% Sofia University "St. Kliment Ohridski"
%
%
% [Z, Z_TYPES, Z_AGES]=simulate_BP(SIM_NUM, T, H, DRAW_S, DRAW_H, DRAW_U, DRAW_Z_0, DRAW_MU, DRAW_IM, APPROX_LIMIT) simulates a
% branching process and returns the number of all particles - Z and the number of particles by type - Z_TYPES.
% Z is a matrix of size NUM_SIM x T/H + 1, containing the simulation paths (the total number of particles for each
% moment of time) in its rows. Z_TYPES is a 3D matrix of size NUM_SIM x N_TYPES x T/H + 1, where N_TYPES is the number
% of types in the branching process. It contains the simulation paths differentiated by types. The sum of all types is
% equal to the total number of particles, stored in Z. Z_AGES contains the simulated age structure of the population 
% by type for every moment of time. It is a 4D matrix of size SIM_NUM x N_AGES x N_TYPES x (T/H+1), where N_AGES is the 
% number of discrete ages, used in the simulation. 
% Note: As a 4D matrix, it might require a lot of RAM. If the age structure is not requested as an output, then it is
% not stored, to save RAM.
%
% SIM_NUM is the number of simulations to perform.
%
% T is the horison of the simulation.
%
% H is the time step for continuous-time branching processes.
%
% DRAW_S could be a matrix sized N_AGES x N_TYPES, representing time-invariant survivability functions for all types of
% particles, where N_AGES is the number of discrete ages in 0:H:<the largest possible age>, and N_TYPES is the number of
% particle types. Or it could be a function that generates the survivability functions for all types of particles on
% every time step. If DRAW_S is a function, then on every call it should generate a random or a constant 3D matrix sized
% N_AGES x N_TYPES x (T/H+1), representing the changing survivability function through time.
%
% Note: If DRAW_S is defined as matrix, the environment is time-invariant.
% If it is defined as a function generating the same N_AGES x N_TYPES x (T/H+1) matrix, then the environment is varying.
% If it generates different N_AGES x N_TYPES x (T/H+1) matrix on every call, then the environment is random.
%
% DRAW_H could be a 2D matrix sized MAX_OFFSPRING x N_TYPES, containing probabilities for 0, 1, ..., MAX_OFFSPRING number
% of children when birth happens, for the different types of particles. It represents the coeficients of the offspring
% generating functions in its columns. If DRAW_H is defined as function, then it should return a 3D matrix sized
% MAX_OFFSPRING x N_TYPES x (T/H+1) on every call. As for DRAW_S, if the matrix is constant, then the environment is
% varying, if it is random, the environment is random.
%
% DRAW_U could be a 2D matrix sized N_TYPES x N_TYPES, containing probabilities for mutations between types. The element
% (i, j) contains the probability for mutation from type j to type i. If DRAW_U is defined as a function, then it should
% return a 3D matrix sized N_TYPES x N_TYPES x (T/H+1) on every call. As for DRAW_S, if the matrix is constant, then the
% environment is varying, if it is random, the environment is random.
%
% DRAW_Z_0 could be a row vector with length N_TYPES, containing the number of initial particles for each type. In that
% case all initial particles are considered to be of age 0 at time 0. Or it could a matrix sized N_AGES x N_TYPES,
% containing the number of initial particles by age. Or it could be a function that generates a 3D matrix sized
% 1 x N_TYPES x (T/H+1) or N_AGES x N_TYPES x (T/H+1), which could be constant or random on every call.
%
% DRAW_MU (Optional) could be a 2D matrix sized N_AGES x N_TYPES, containing the expected number of births a woman has until age x.
% It could also be defined as a function that returns a matrix sized N_AGES x N_TYPES x (T/H+1), which could be constant
% or random on every call. Use DRAW_MU=[] to model Galton-Watson and Bellman-Harris branching processes, for which
% no point process is defined. 
%
% DRAW_IM (Optional) could be a 2D matrix sized N_AGES x N_TYPES, containing the number of immigrants on all ages by type.
% It could also be defined as a function that returns a matrix sized N_AGES x N_TYPES x (T/H+1), which could be constant
% or random on every call. When DRAW_IM=[], there is no immigration.
%
% APPROX_LIMIT (Optional, Default=20) determines when the normal approximation of binomial
% distribution is appropriate. When N*p>=approx_limit && N*(1-p)>=approx_limit, the algorithm uses normal approximation.
% This allows handling of large populations without consuming more RAM or CPU time.
% 
% VIRTUES AND RESTRICTIONS OF THE SIMULATION METHOD
%   Theoretical Virtues:
%    - It is one single code that simulates Galton-Watson, Bellman-Harris and General Branching Processes;
%    - It cnsiders multi-type processes with mutations between them;
%    - Mutation probabilities could be constant, varying or random;
%    - It considers constant, varying or random environment;
%    - It considers constant, varying or random immigration;
%    - It considers constant, varying or random initial number of particles on constant, varying or random ages;
%    - It produces not only the total population count, but also the population count by age and type at each moment of 
%       time. 
%    - It could be extended to include controlled branching processes. However, this would most probably require the 
%       code to be customized only for the specific process as the possibilities for such theoretical dependences 
%       between population size, birth, death and migration laws could be quite large;
%   
%   Computational Virtues:
%    - It presents a simple and short code that simulates branching processes;
%    - It is capable of simulating VERY large number of particles (for example: 10^250) in the branching processes, 
%       without requiring a lot of RAM (by using normal approximation of binomial distribution when possible);
%    - The simulation is actually faster for larger populations due to the normal approximation of multinomial 
%       distribution;
%    - It uses all computer cores for faster calculation.
%
%   Restrictions:
%    - Particles cannot die/give birth at the beginning of their life. Probability of that event is considered zero;
%    - Birth and death densities must be smooth functions with exception of finite jump type discontinuities;
%    - The birth, death distributions, mutation probabilities and immigration are independent with each other and with 
%       the branching process itself. So the simulation does not include the class of "controlled branching processes",
%       although it could be extended to suit the specific needs. Depending on the type of controlled process and the 
%       dependence on the age structure at previous times, this task could require a lot of RAM;
%    - It considers only immigration, as emigration is naturally dependent on the population count and age structure 
%       and could be defined in a variety of ways.
%    - Returning the population count by age as an output may require a lot of computer memory (more than a 100 GB RAM
%       in some cases);
%
%
%
%       BranchingProcessSimulator - opensource program for simulating branching processes
%       Copyright (C) 2018  Plamen Ivaylov Trayanov
%
%       This program is free software: you can redistribute it and/or modify
%       it under the terms of the GNU General Public License as published by
%       the Free Software Foundation, either version 3 of the License, or
%       (at your option) any later version.
% 
%       This program is distributed in the hope that it will be useful,
%       but WITHOUT ANY WARRANTY; without even the implied warranty of
%       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%       GNU General Public License for more details.
% 
%       You should have received a copy of the GNU General Public License
%       along with this program. If not, see <http://www.gnu.org/licenses/>.
% 


%% input parser, basic check and adjustment of the input format
Results=customInputParser(sim_num, T, h, draw_S, draw_H, draw_U, draw_Z_0, varargin{:});

getAgeStructure=(nargout==3);   % determines whether the age structure is requested or not
T=ceil(T/h)*h;      % make T divisible by h

%% perform the simulation using all cores
n_types=size(Results.draw_S(),2);    % the number of types
Z_types=zeros(sim_num, n_types, T/h+1);     % the branching process by types - sum of all individuals by types at time [0, T].

if getAgeStructure      % setting it to false, saves RAM, as Z_age requires a lot of memory
    Z_ages=zeros([sim_num, size(Results.draw_S())]);     % population count by [sim_num, age, type, time]
end
parfor ind=1:sim_num
    % for every simulation path, before simulating the branching process, generate simulation paths of the probability 
    % laws that determine the branching process properties. Later they are used for the branching process simulation.
    S=Results.draw_S();
    H=Results.draw_H();
    U=Results.draw_U();
    Z_0=Results.draw_Z_0();
    if ~isempty(Results.draw_mu)    % draw mu if it is specified in the input, i.e. the process is general branching process
        mu=Results.draw_mu(); 
    end
    if ~isempty(Results.draw_Im)    % draw Im if it is specified in the input
        Im=Results.draw_Im(); 
    end
    
    % initialize the array in which we store the age structure of the
    % population. Instead of saving it for every simulation and returning it as output, here we
    % initialize it. This is done to reduce the amount of RAM needed to
    % store the output
    N=zeros(size(S));     % population structure - N(age, type, time)
    if(size(Z_0,1)~=size(S,1))  % if Z_0 does not contain information for ages
        N(1,:,1)=Z_0;   % then assume the given initial population is of age 0
    else
        N(:,:,1)=Z_0;   % then load the given initial population age structure
    end
    
    for t=1:T/h       % go through time
        if ~isempty(Results.draw_Im)  % add immigration at the beginning of the time interval 
            % (this also allows the initial population to be defined as immigrants at time 0)
            N(:,:,t)=N(:,:,t)+Im(:,:,t);
        end
        % if there is at least one particle of any type at time t, simulate for t+h
        if any(any(N(:,:,t))) || ~isempty(Results.draw_Im)
            for i=1:size(N,1)-1    % for each age
                for j=1:size(N,2)  % for each type of particle
                    if N(i,j,t)~=0 && S(i,j,t)~=0   % if no particles of that type are alive, no need to simulate their deaths
                        % simulate the number of deaths at time t
                        deaths=N(i,j,t)-binornd_large(N(i,j,t), S(i+1,j,t)/S(i,j,t), Results.approx_limit);
                        if isempty(Results.draw_mu)     % assumes the particle splits at its death if mu is unspecified
                            % BGW, BH branching process - the births happen when death occures
                            births=mnrnd_large(deaths, H(:,j,t)',1, Results.approx_limit)*(0:(length(H(:,j,t))-1))';
                        else
                            % the case of General Branching Process - the births occur according to a point process,
                            % during the life of the particles
                            births=mnrnd_large(binornd_large(N(i,j,t), mu(i,j,t)*h, Results.approx_limit), H(:,j,t)',1, Results.approx_limit)*(0:(length(H(:,j,t))-1))';
                        end
                        N(i+1,j,t+1)=N(i,j,t)-deaths;   % the ones that survived get to get older by h
                        N(1,:,t+1)=N(1,:,t+1) + mnrnd_large(births, U(:,j,t)',1, Results.approx_limit);     % simulate the mutations to other types
                    end
                end
            end
        else
            % if the population count is 0 at time t, and we have no immigration, then no need to simulate any further
            break;
        end
    end
    if getAgeStructure  % if we want the age structure, save it in Z_ages
        Z_ages(ind, :, :, :)=N;
    end
    Z_types(ind, :, :)=sum(N,1);  % sum of all individuals by type (no age information)
end
Z=squeeze(sum(Z_types, 2));     % sum of all individuals, no age or type information
end

%% private functions used by the simulator
function Results=customInputParser(sim_num, T, h, draw_S, draw_H, draw_U, draw_Z_0, varargin)
p=inputParser();
p.addRequired('sim_num', @(x)(isnumeric(x) && isscalar(x) && x>0 && mod(x,1)==0)); % check if sim_num is a positive integer number
p.addRequired('T', @(x)(isnumeric(x) && isscalar(x) && x>0)); % check if T is a positive number
p.addRequired('h', @(x)(isnumeric(x) && isscalar(x) && x>0)); % check if h is a positive number

% check if draw_S is a numerical matrix of decreasing numbers in the columns, between 0 and 1
% OR a function that returns a sequence of such matricies.
validateS_matrix=@(x)(isnumeric(x) && ismatrix(x) && all(x(:)>=0 & x(:)<=1) && all(all(diff(x)<=0)));
p.addRequired('draw_S', @(x)(validateS_matrix(x) || (isa(x,'function_handle') && isnumeric(x()))));
validate_ProbabilityMatrix=@(x)(isnumeric(x) && all(x(:)>=0 & x(:)<=1));
p.addRequired('draw_H', @(x)(validate_ProbabilityMatrix(x) || (isa(x,'function_handle') && isnumeric(x()))));
p.addRequired('draw_U', @(x)(validate_ProbabilityMatrix(x) || (isa(x,'function_handle') && isnumeric(x()))));
p.addRequired('draw_Z_0', @(x)((isnumeric(x) && all(x(:)>=0 & mod(x(:),1)==0)) || (isa(x,'function_handle') && isnumeric(x()))));
p.addOptional('draw_mu', [], @(x)((isnumeric(x) && ismatrix(x) && all(x(:)>=0)) || ...
    (isa(x,'function_handle') && isnumeric(x()))));
p.addOptional('draw_Im', [], @(x)((isnumeric(x) && ismatrix(x) && all(x(:)>=0 & mod(x(:),1)==0)) || ...
    (isa(x,'function_handle') && isnumeric(x()))));
p.addOptional('approx_limit', 20, @(x)(isnumeric(x) && isscalar(x) && x>0)); % check if T is a positive number

p.parse(sim_num, T, h, draw_S, draw_H, draw_U, draw_Z_0, varargin{:});
Results = p.Results;

%% checks the format of the input. If the input is in matrix form, transform it to a function form.
if(isnumeric(draw_S))
    if ndims(draw_S) == 3   % the case of varying (but not random) environment, i.e. draw_S contains survivability 
                            % function for every moment of time
        Results.draw_S=@()(draw_S);     % use different but non-random survivability functions
    else    % the case of constant environment 
        % when the survivability function is a 2D matrix, make the matrix the same for all moments of time 0:h:T,
        % i.e. the survivability function is time invariant
        Results.draw_S=@()(repmat(draw_S, [1,1,T/h+1]));
    end
end
if(isnumeric(draw_H))   % do the same as for draw_S
    if ndims(draw_H) == 3
        Results.draw_H=@()(draw_H);
    else
        Results.draw_H=@()(repmat(draw_H,[1,1,T/h+1]));
    end
end
if(isnumeric(draw_U))   % do the same as for draw_S
    if ndims(draw_U) == 3
        Results.draw_U=@()(draw_U);
    else
        Results.draw_U=@()(repmat(draw_U,[1,1,T/h+1]));
    end
end
if(~isempty(Results.draw_mu) && isnumeric(Results.draw_mu))   % do the same as for draw_S
    if ndims(Results.draw_mu) == 3
        Results.draw_mu=@()(Results.draw_mu);
    else
        Results.draw_mu=@()(repmat(Results.draw_mu,[1,1,T/h+1]));
    end
end

if(~isempty(Results.draw_Im) && isnumeric(Results.draw_Im))   % do the same as for draw_S
    if ndims(Results.draw_Im) == 3
        Results.draw_Im=@()(Results.draw_Im);
    else
        Results.draw_Im=@()(repmat(Results.draw_Im,[1,1,T/h+1]));
    end
end

if(isnumeric(Results.draw_Z_0))
    Results.draw_Z_0=@()(Results.draw_Z_0);
end
end

function res=binornd_large(N, p, approx_limit)
% simulation from the binomial distribution for large numbers
    if p==0     % if the probability is 0, then no need to simulate (saves time)
        res=0;
    elseif p==1     % if p is 1, then again, no need to simulate
        res=N;
    elseif N*p<approx_limit || N*(1-p)<approx_limit     % criteria for using normal approximation
        res=binornd(N,p);   % use binomial distribution
    else
        % use the normal approximation
        res=round(normrnd(N*p, sqrt(N*p*(1-p))));
    end
end

function res=mnrnd_large(N, p, v, approx_limit)
    % simulate from the multinomial distribution for small and large samples
    res=zeros(size(p));
    i_pos=p~=0;     % indecies of the positive probabilities
    p_pos=p(i_pos);         % The positive propabilities
    if length(p_pos)==1     % only 1 positive probability, i.e outcome is certain, so do not simulate
        res(i_pos)=N;       % the only positive probability should be 1
    elseif any(N*p_pos<approx_limit | N*(1-p_pos)<approx_limit)     % criteria for using normal approximation
        res(i_pos)=mnrnd(N, p_pos, v);      % get multinomial simulations only for the positive probabilities
    else
        % use the normal approximation, i.e. multivariate normal distribution with specific covariance matrix,
        % which makes the sum of the simulated values equal to N
        Sigma=diag(p_pos)-p_pos'*p_pos;
        res(i_pos)=round(sqrt(N)*mvnrnd(sqrt(N)*p_pos, Sigma, v));
    end
end
