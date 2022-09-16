%Analyze all networks for the vary off add k 700 parameter ensemble using
%diferent thresholds for binarizing the data to check for robustness

addpath(genpath(pwd))
addpath(genpath('./../Functions'))
addpath(genpath('./../Data'))

clear
clc

save_directory = './../../../../Paper/extractedData/20210528_for_paper/binarization_threshold_test/';
n_params = 700;
trim_length = 100;
n_alleles = 2;


for irep = 1
    irep
    runID = sprintf('20210604_vary_off_add_k_for_paper_long-rep%d', irep);

    min_species = 5;
    max_species = 5;
    export_bursts=[];
    export_OR=[];
    export_upstream_bursts=[];

    for ispecies = min_species:max_species
        ispecies
        load(sprintf('M_iso%d', ispecies))
        max_subnet = size(M_iso,2);

        for isubnet = 1:max_subnet

            load(sprintf('S_outpar_%s_%d_%d_1', runID, ispecies, isubnet));
            
            for threshold = 1:5

                for iparam = 1:n_params

                    [bursts, bursts_duration, bursts_total_nodes, bursts_concurrent_nodes, nodes_fraction_on]=...
                        find_RNA_bursts(ispecies,n_alleles, threshold, S_outpar{iparam}, trim_length, 1);

                    [mean_odds_ratio, all_odds_ratio]=...
                        find_RNA_allele_OR(ispecies,n_alleles,threshold,S_outpar{iparam},trim_length);

                    export_bursts = [export_bursts; ispecies isubnet iparam irep threshold length(bursts) mean(bursts_duration) mean(bursts_total_nodes) mean(bursts_concurrent_nodes) mean(mean(nodes_fraction_on,1))];
                    export_OR = [export_OR; ispecies isubnet iparam irep threshold mean_odds_ratio];
                end
            end
        end
    end
    writematrix(export_bursts, strcat(save_directory, sprintf('long_bursts-rep%d.csv', irep)));
    writematrix(export_OR, strcat(save_directory, sprintf('long_OR-rep%d.csv', irep)));
end

load('20210407_vary_off_add_k_fewer_params_for_signalburst')
writematrix(DataParams, strcat(save_directory, 'metadata.csv'));