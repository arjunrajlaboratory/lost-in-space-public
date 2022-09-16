%Analyze the 2 node, 1 way custom network, A -> B
% Ultimately not used in Paper

addpath(genpath(pwd))
addpath(genpath('./../Functions'))
addpath(genpath('./../Data'))

clear
clc

save_directory = './../../../../Paper/extractedData/20210528_for_paper/2node_1way/';

runID = '20210527_vary_off_add_k_for_paper_1way_2node-rep1';
n_species = 2;
n_alleles = 2;
threshold = 3;
n_params = 700;
trim_length = 100;

load(sprintf('S_outpar_%s_2_1_1', runID));
export_bursts=[];
export_OR=[];
export_upstream_bursts=[];
for iparam = 1:n_params
    
    [bursts, bursts_duration, bursts_total_nodes, bursts_concurrent_nodes, nodes_fraction_on]=...
        find_RNA_bursts(n_species,n_alleles, threshold, S_outpar{iparam}, trim_length, 1);
    
    %Calc bursts of upstream factor for normalization
    [upstream_bursts, upstream_bursts_duration] = ...
        find_upstream_bursts(n_species,n_alleles, threshold, S_outpar{iparam}, trim_length);

    
    [mean_odds_ratio, all_odds_ratio]=...
        find_RNA_allele_OR(n_species,n_alleles,threshold,S_outpar{iparam},trim_length);
    
    export_bursts = [export_bursts; n_species iparam length(bursts) mean(bursts_duration) mean(bursts_total_nodes) mean(bursts_concurrent_nodes) mean(nodes_fraction_on,1)];
    export_OR = [export_OR; n_species iparam mean_odds_ratio all_odds_ratio'];
    export_upstream_bursts = [export_upstream_bursts; n_species iparam length(upstream_bursts) mean(upstream_bursts_duration)];
end


writematrix(export_bursts, strcat(save_directory, 'bursts.csv'));
writematrix(export_upstream_bursts, strcat(save_directory, 'upstream_bursts.csv'));
writematrix(export_OR, strcat(save_directory, 'OR.csv'));

load('20210407_vary_off_add_k_fewer_params_for_signalburst')
writematrix(DataParams, strcat(save_directory, 'metadata.csv'));