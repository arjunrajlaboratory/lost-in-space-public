% Generate data the 2-5 node network, all subnets, vary off add k 700 parameter ensemble
% with the longer run length used in most analyses in the paper
% Run this first since it also generates the 700 parameter ensemble used in
% most other simulations

addpath(genpath(pwd))
addpath(genpath('./../Functions'))
addpath(genpath('./../Data'))
clear
clc

% Generate 700 parameter ensemble
r_deg1_to_test = .01;
r_prod1_ratios = .01;
r_off1_ratios = logspace(0, 2, 10);
r_add1_ratios = logspace(-1, 2, 10);
r_on1_ratios = .025;

proddiff1 = 10000;
k1 = linspace(20, 200, 7);

n1 = 1;
p = 1; % scaling factor not ultimately used in final version

combos = combvec(r_prod1_ratios,...
    r_deg1_to_test, r_add1_ratios, ...
    r_off1_ratios, proddiff1,...
    r_on1_ratios, k1, n1, p)';

corrected_combos = [combos(:,1).*combos(:,2) combos(:,2) combos(:,3).*combos(:,2).*combos(:,9)...
    combos(:,4).*combos(:,2).*combos(:,9) combos(:,5) combos(:,6).*combos(:,2).*combos(:,9) combos(:,7:8) combos(:,9)];

DataParams = corrected_combos;

save('./../Data/20210407_vary_off_add_k_fewer_params_for_signalburst','DataParams')


% run simulation on 2-5 node networks. Networks larger than 5 nodes were
% too big to run with this simulation length

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
    runID = sprintf('20210604_vary_off_add_k_for_paper_long-rep%d', 1); % unique identifier shared by all output files
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    generateData_allele_diffRNA(nruns, n_species, n_alleles, maxgillespie/100, maxgillespie,gen,init,data, type, rand_gen, rand_set, runID);
end
