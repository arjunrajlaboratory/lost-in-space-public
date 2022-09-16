% Generate data the 5 node network, vary off add n 500 parameter ensemble
% with the longer run length used in most analyses in the paper
% This data set is only used to verify robustness to changes in n and so is
% only run with the 5 node network

addpath(genpath(pwd))
addpath(genpath('./../Functions'))
addpath(genpath('./../Data'))

clear
clc

% Generate 500 parameter ensemble
r_deg1_to_test = .01;
r_prod1_ratios = .01;
r_off1_ratios = logspace(0, 2, 10);
r_add1_ratios = logspace(-1, 2, 10);
r_on1_ratios = .025;


proddiff1 = 10000;
k1 = 110;

n1 = [0.2 0.5 1 2 4];

p = 1;  %scaling factor not ultimately used in final version
combos = combvec(r_prod1_ratios,...
    r_deg1_to_test, r_add1_ratios, ...
    r_off1_ratios, proddiff1,...
    r_on1_ratios, k1, n1, p)';

corrected_combos = [combos(:,1).*combos(:,2) combos(:,2) combos(:,3).*combos(:,2).*combos(:,9)...
    combos(:,4).*combos(:,2).*combos(:,9) combos(:,5) combos(:,6).*combos(:,2).*combos(:,9) combos(:,7:8) combos(:,9)];
DataParams = corrected_combos;
save('./../Data/20211112_vary_off_add_n_params','DataParams')


clear
clc
% run simulation on 5 node networks. 
for ispecies = 5
    for irep = 1
    nruns = 500;                             %corresponds to number of parameter sets to be run
    n_species = ispecies;                           %number of nodes
    n_alleles = 2;
    maxgillespie = 1000000*20;                 %number of Gillespie simulated time units
    gen = 'no';                             %yes = generate new parameters (or 'no' = load parameters)
    init.spec = repmat(0,1,(n_species*n_alleles));     %initial values for all species (at t=0)
    init.Bon = zeros(1,(n_species*n_alleles));          %initial values for the burst (where 0 = 'off' and 1 = 'on')
    data = '20211112_vary_off_add_n_params';
    type = 'normal';
    rand_gen = 'no';
    rand_set = irep;
    runID = sprintf('20211112_vary_off_add_n_params_long-rep%d', irep);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    generateData_allele_diffRNA(nruns, n_species, n_alleles, maxgillespie/100, maxgillespie,gen,init,data, type, rand_gen, rand_set, runID);
    end
end
