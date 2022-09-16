%This script generates the simple 2 node asymmetric network given in figure 1
%and generates data to be plotted in figure 1

addpath(genpath(pwd))
addpath(genpath('./../Functions'))
addpath(genpath('./../Data'))
clear
clc

%Make custom 1 way 2 node network
load('M_iso2')
inet = M_iso{1};
inet(2,1) = 0;
M_iso = {inet};
save('./../Data/M_iso2_custom', 'M_iso');

clear
clc

%Generate data
for irep = 1
    nruns = 700;                             %corresponds to number of parameter sets to be run
    n_species = 2;                           %number of nodes
    n_alleles = 2;                          %number of alleles
    maxgillespie = 1000000;                 %number of Gillespie simulated time units
    gen = 'no';                             %yes = generate new parameters (or 'no' = load parameters)
    init.spec = repmat(0,1,(n_species*n_alleles));     %initial values for all species (at t=0)
    init.Bon = zeros(1,(n_species*n_alleles));          %initial values for the burst (where 0 = 'off' and 1 = 'on')
    data = '20210407_vary_off_add_k_fewer_params_for_signalburst'; % Parameters to load
    type = 'normal';                                            % 'normal' for symmetric parameters
    rand_gen = 'no';                                             % 'yes' to generate a random seed. 'no' to use rand_set for seed
    rand_set = irep;                                            % RNG seed, only used if rand_gen = 'no'
    runID = sprintf('20210527_vary_off_add_k_for_paper_1way_2node-rep%d', irep); % unique identifier shared by all output files
    
    generateData_allele_diffRNA_custom_1way_2node(nruns, n_species, n_alleles, maxgillespie/100, maxgillespie,gen,init,data, type, rand_gen, rand_set, runID);
end