% Generate data the linear 2-5 node network, vary off add k 700 parameter ensemble
% with the longer run length used in most analyses in the paper

addpath(genpath(pwd))
addpath(genpath('./../Functions'))
addpath(genpath('./../Data'))
clear
clc


for ispecies = 2:5
    nruns = 700;                             %corresponds to number of parameter sets to be run
    n_species = ispecies;                           %number of nodes
    n_alleles = 2;                           % number of alleles
    maxgillespie = 1000000*20;                 %number of Gillespie simulated time units
    gen = 'no';                             %yes = generate new parameters (or 'no' = load parameters)
    init.spec = repmat(0,1,(n_species*n_alleles));     %initial values for all species (at t=0)
    init.Bon = zeros(1,(n_species*n_alleles));          %initial values for the burst (where 0 = 'off' and 1 = 'on')
    data = '20210407_vary_off_add_k_fewer_params_for_signalburst'; % Parameters to load
    type = 'normal';                                    % 'normal' for symmetric parameters
    rand_gen = 'no';                                    % 'yes' to generate a random seed. 'no' to use rand_set for seed
    rand_set = 1;                                       % RNG seed, only used if rand_gen = 'no
    runID = sprintf('revision_vary_off_add_k_for_paper_long_linear_network-rep%d', 1); % unique identifier shared by all output files
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    generateData_allele_diffRNA_linear(nruns, n_species, n_alleles, maxgillespie/100, maxgillespie,gen,init,data, type, rand_gen, rand_set, runID);
end
