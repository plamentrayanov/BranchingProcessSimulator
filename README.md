# BranchingProcessSimulator
Simulates multi-type Galton-Watson, Bellman-Harris and Crump-Mode-Jagers branching processes with immigration - in constant, varying or random environment. The process is allowed to start from a random number of particles on different ages. The mutation probabilities are also allowed to be random. The features and restrictions of the simulation method are described below:

FEATURES
Theoretical features:
- One single code that simulates Galton-Watson, Bellman-Harris and General Branching Processes;
- It considers multi-type processes with mutations between them;
- Mutation probabilities could be constant, varying or random;
- It considers constant, varying or random environment;
- It considers constant, varying or random immigration;
- It considers constant, varying or random initial number of particles on constant, varying or random ages;
- It produces the total population count and the population count by age and type at each moment of time.
- It could be extended to include controlled branching processes. However, this would most probably require the code to be customized only for the specific process as the possibilities for such theoretical dependences between population size, birth, death and migration laws could be quite large;

Computational features:
- It presents a simple and short code that simulates branching processes;
- It is capable of simulating VERY large number of particles (for example: 10^250) in the branching processes, without requiring a lot of RAM (by using normal approximation of binomial distribution when possible);
- The simulation is actually faster for larger populations due to the normal approximation of multinomial distribution;
- It uses all computer cores for faster calculation.

Restrictions:
- Particles cannot die/give birth at the beginning of their life. Probability of that event is considered zero;
- Birth and death densities must be smooth functions with exception of finite jump type discontinuities;
- The birth, death distributions, mutation probabilities and immigration are independent with each other and with the branching process itself. So the simulation does not include the class of "controlled branching processes", although it could be extended to suit the specific needs. Depending on the type of controlled process and the dependence on the age structure at previous times, this task could require a lot of RAM;
- It considers only immigration, as emigration is naturally dependent on the population count and age structure and could be defined in a variety of ways.
- Returning the population count by age as an output may require a lot of computer memory (more than a 100 GB RAM in some cases);
