% Generate data the 2-5 node network, all subnets, vary off add k 700 parameter ensemble
% with a shorter run length for instances where a smaller individual file
% size is useful for plotting. Generate 15 replicates (used in a
% since replaced analysis that pooled replicates)

addpath(genpath(pwd))
addpath(genpath('./../Functions'))
addpath(genpath('./../Data'))
clear
clc

for ispecies = 2:5
    for irep = 1:15
    nruns = 700;                             %corresponds to number of parameter sets to be run
    n_species = ispecies;                           %number of nodes
    n_alleles = 2;
    maxgillespie = 1000000;                 %number of Gillespie simulated time units
    gen = 'no';                             %yes = generate new parameters (or 'no' = load parameters)
    init.spec = repmat(0,1,(n_species*n_alleles));     %initial values for all species (at t=0)
    init.Bon = zeros(1,(n_species*n_alleles));          %initial values for the burst (where 0 = 'off' and 1 = 'on')
    data = '20210407_vary_off_add_k_fewer_params_for_signalburst'; % Parameters to load
    type = 'normal';                                    % 'normal' for symmetric parameters
    rand_gen = 'no';                                    % 'yes' to generate a random seed. 'no' to use rand_set for seed
    rand_set = irep;                                    % RNG seed, only used if rand_gen = 'no'
    runID = sprintf('20210524_vary_off_add_k_for_paper-rep%d', irep);       % unique identifier shared by all output files
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    generateData_allele_diffRNA(nruns, n_species, n_alleles, maxgillespie/100, maxgillespie,gen,init,data, type, rand_gen, rand_set, runID);
    end
end
