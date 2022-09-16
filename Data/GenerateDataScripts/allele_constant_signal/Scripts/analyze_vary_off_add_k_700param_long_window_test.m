%Find bursts using different 'burst windows' to check for robustness in
%response time v burst metrics 

addpath(genpath(pwd))
addpath(genpath('./../Functions'))
addpath(genpath('./../Data'))

clear
clc

save_directory = './../../../../Paper/extractedData/20210528_for_paper/burst_window_test/';
n_params = 700;
trim_length = 100;
n_alleles = 2;
threshold=3;

for irep = 1
    irep
    runID = sprintf('20210604_vary_off_add_k_for_paper_long-rep%d', irep);

    min_species = 5;
    max_species = 5;
    export_bursts=[];

    for ispecies = min_species:max_species
        ispecies
        load(sprintf('M_iso%d', ispecies))
        max_subnet=1; %only use subnet=1 since that is the subnetwork on which signal responsiveness was simulated
        for isubnet = 1:max_subnet

            load(sprintf('S_outpar_%s_%d_%d_1', runID, ispecies, isubnet));
            
            for burst_window = 2:6

                for iparam = 1:n_params

                    [bursts, bursts_duration, bursts_total_nodes, bursts_concurrent_nodes, nodes_fraction_on]=...
                        find_RNA_bursts_with_window(ispecies,n_alleles, threshold, S_outpar{iparam}, trim_length, 1, burst_window);

                    export_bursts = [export_bursts; ispecies isubnet iparam irep burst_window length(bursts) mean(bursts_duration) mean(bursts_total_nodes) mean(bursts_concurrent_nodes) mean(mean(nodes_fraction_on,1))];
                end
            end
        end
    end
    writematrix(export_bursts, strcat(save_directory, sprintf('long_bursts-rep%d.csv', irep)));
end

load('20210407_vary_off_add_k_fewer_params_for_signalburst')
writematrix(DataParams, strcat(save_directory, 'metadata.csv'));