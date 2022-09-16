%Classify negative regulation networks as constitutively expressing or not in order to
%filer in downstream plotting and analysis

addpath(genpath(pwd))
addpath(genpath('./../Functions'))
addpath(genpath('./../Data'))

clear
clc
save_directory = './../../../../Paper/extractedData/20210528_for_paper/negative_regulation_networks/';
threshold = 3; %threshold with which to binarize data
n_params = 700;
trim_length = 100; %length to trim off the start of simulation to account for initialization time
n_alleles = 2;

upper_bound = 0.9995; %fraction of time simulation needs to be above threshold to call as constitutively expressing

for irep = 1
    irep
    runID = sprintf('revision_negative_regulation_for_paper_long-rep%d', irep);
    
    min_species = 5;
    max_species = 5;
    export_sim_class=[];
    
    for ispecies = min_species:max_species
        ispecies
        max_subnet = 3;
        
        for isubnet = 1:max_subnet
        
        load(sprintf('S_outpar_%s_%d_%d_1', runID, ispecies, isubnet));
        
        
        for iparam = 1:n_params
            [fraction_all_above, is_constant_high] = find_fraction_all_above_threshold(ispecies,n_alleles, threshold, S_outpar{iparam}, trim_length, upper_bound);
            export_sim_class = [export_sim_class; ispecies isubnet iparam irep fraction_all_above,is_constant_high];
        end
        end
    end
    writematrix(export_sim_class, strcat(save_directory, sprintf('long_simulation_class-rep%d.csv', irep)));
end