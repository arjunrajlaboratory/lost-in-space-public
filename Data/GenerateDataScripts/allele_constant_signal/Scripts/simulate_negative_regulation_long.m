% Generate data for custom 5 node networks with negative regulation using
% custom 700 parameter ensemble with the longer run length used in most analyses in the paper

addpath(genpath(pwd))
addpath(genpath('./../Functions'))
addpath(genpath('./../Data'))

clear
clc

for ispecies = 5
    nruns = 700;                             %corresponds to number of parameter sets to be run
    n_species = ispecies;                           %number of nodes
    n_alleles = 2;                           % number of alleles
    maxgillespie = 1000000*20;                 %number of Gillespie simulated time units
    gen = 'no';                             %yes = generate new parameters (or 'no' = load parameters)
    init.spec = repmat(0,1,(n_species*n_alleles));     %initial values for all species (at t=0)
    init.Bon = zeros(1,(n_species*n_alleles));          %initial values for the burst (where 0 = 'off' and 1 = 'on')
    data = '20220318_negative_regulation_params'; % Parameters to load
    type = 'normal';                                    % 'normal' for symmetric parameters
    rand_gen = 'no';                                    % 'yes' to generate a random seed. 'no' to use rand_set for seed
    rand_set = 1;                                       % RNG seed, only used if rand_gen = 'no
    runID = sprintf('revision_negative_regulation_for_paper_long-rep%d', 1); % unique identifier shared by all output files
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %uses custom generateData function to allow for custom negative
    %regulation networks
    generateData_allele_diffRNA_neg_reg(nruns, n_species, n_alleles, maxgillespie/100, maxgillespie,gen,init,data, type, rand_gen, rand_set, runID);
end
