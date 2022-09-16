%Classify networks as constitutively expressing or not in order to
%filter in downstream plotting and analysis

addpath(genpath(pwd))
addpath(genpath('./../Functions'))
addpath(genpath('./../Data'))

clear
clc
save_directory = './../../../../Paper/extractedData/20210528_for_paper/vary_off_add_k_700params/';
threshold = 3;
n_params = 700;
trim_length = 100;
n_alleles = 2;

upper_bound = 0.9995;

for irep = 1
    irep
    runID = sprintf('20210604_vary_off_add_k_for_paper_long-rep%d', irep);

    min_species = 2;
    max_species = 5;
    export_sim_class=[];

    for ispecies = min_species:max_species
        ispecies
        load(sprintf('M_iso%d', ispecies))
        max_subnet = size(M_iso,2);

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