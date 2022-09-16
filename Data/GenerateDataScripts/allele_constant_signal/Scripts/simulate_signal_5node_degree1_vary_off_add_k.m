addpath(genpath(pwd))
addpath(genpath('./../Functions'))
addpath(genpath('./../Data'))
clear
clc

% Constantly high signal 400 molecules

% Prerun 100 times then run 100 replicates from each

signal_degs = 0;
for ipre_rep = 1:100
    
nruns = 700;                             %corresponds to number of parameter sets to be run
n_species = 5;                           %number of nodes
n_alleles = 2;
maxgillespie = 20000;                 %number of Gillespie simulated time units
gen = 'no';                             %yes = generate new parameters (or 'no' = load parameters)
init.spec = repmat(0,1,(n_species*n_alleles));     %initial values for all species (at t=0)
init.Bon = zeros(1,(n_species*n_alleles));          %initial values for the burst (where 0 = 'off' and 1 = 'on')
data = '20210407_vary_off_add_k_fewer_params_for_signalburst';
type = 'normal';
prerunID = sprintf('20210518_vary_off_add_k_for_constant_signal-prerep%d', ipre_rep);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

generate_prerun_allele(nruns, n_species, n_alleles, maxgillespie/100, maxgillespie,gen,init,data, type, prerunID);


    signal_value = 400;
    num_replicates = 100;
    nruns = 700;                             %corresponds to number of parameter sets to be run
    n_species = 5;                           %number of nodes
    n_alleles = 2;
    maxgillespie = 20000;                 %number of Gillespie simulated time units
    gen = 'no';                             %yes = generate new parameters (or 'no' = load parameters)
    init.spec = repmat(0,1,(n_species*n_alleles));     %initial values for all species (at t=0)
    init.Bon = zeros(1,(n_species*n_alleles));          %initial values for the burst (where 0 = 'off' and 1 = 'on')
    data = '20210407_vary_off_add_k_fewer_params_for_signalburst';
    type = 'normal';
    savefileID = sprintf('20210518_vary_off_add_k_for_constant_signal_%d-prerep%d', signal_value, ipre_rep);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    generateData_signal_allele_multiplereps(nruns, n_species, n_alleles, maxgillespie/100, maxgillespie,gen,init,data, type, savefileID, signal_value, signal_degs, prerunID, num_replicates)
end
