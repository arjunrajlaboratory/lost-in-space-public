%Revision script
%Analyze all linear networks using original metrics and revision
%cross-correlation metrics

addpath(genpath(pwd))
addpath(genpath('./../Functions'))
addpath(genpath('./../Data'))

clear
clc
save_directory = './../../../../Paper/extractedData/20210528_for_paper/linear_networks/';
threshold = 3; %threshold with which to binarize data
n_params = 700;
trim_length = 100; %length to trim off the start of simulation to account for initialization time
n_alleles = 2;


for irep = 1 %the 'long' simulations only use one replicate
    irep
    runID = sprintf('revision_vary_off_add_k_for_paper_long_linear_network-rep%d', irep);
    
    min_species = 2;
    max_species = 8;
    
    
    for ispecies = min_species:max_species
        ispecies
        export_bursts=[];
        export_OR=[];
        isubnet = 1;
        export_rows = ispecies * (ispecies - 1) / 2;
        all_xcorrs_bursts_export=[];
        max_xcorrs_bursts_export=[];
        all_xcorrs_allele_export=[];
        max_xcorrs_allele_export=[];
        export_per_gene_above=[];
        export_bursts=[];
        load(sprintf('S_outpar_%s_%d_%d_1', runID, ispecies, isubnet));
        
        
        for iparam = 1:n_params
            
            [per_gene_above_threshold]=...
                find_above_threshold_per_gene(ispecies,n_alleles, threshold, S_outpar{iparam}, trim_length);
            
            [all_xcorrs_bursts, max_xcorrs_bursts]=...
                find_xcorrs_RNA_bursts(ispecies,n_alleles,S_outpar{iparam},trim_length);

            [all_xcorrs_allele, max_xcorrs_allele]=...
                find_xcorrs_alleles(ispecies,n_alleles,S_outpar{iparam},trim_length);

            [mean_odds_ratio, all_odds_ratio]=...
                find_RNA_allele_OR(ispecies,n_alleles,threshold,S_outpar{iparam},trim_length);
            
            [bursts, bursts_duration, bursts_total_nodes, bursts_concurrent_nodes, nodes_fraction_on]=...
                find_RNA_bursts(ispecies,n_alleles, threshold, S_outpar{iparam}, trim_length, 1);
            
            all_xcorrs_bursts_export = [all_xcorrs_bursts_export;...
                repelem(ispecies, export_rows)'...
                repelem(isubnet, export_rows)'...
                repelem(iparam, export_rows)'...
                repelem(irep, export_rows)'...
                all_xcorrs_bursts];

            max_xcorrs_bursts_export =[max_xcorrs_bursts_export;...
                repelem(ispecies, export_rows)'...
                repelem(isubnet, export_rows)'...
                repelem(iparam, export_rows)'...
                repelem(irep, export_rows)'...
                max_xcorrs_bursts];

            all_xcorrs_allele_export = [all_xcorrs_allele_export;...
                repelem(ispecies, ispecies)'...
                repelem(isubnet, ispecies)'...
                repelem(iparam, ispecies)'...
                repelem(irep, ispecies)'...
                all_xcorrs_allele];

            max_xcorrs_allele_export =[max_xcorrs_allele_export;...
                repelem(ispecies, ispecies)'...
                repelem(isubnet, ispecies)'...
                repelem(iparam, ispecies)'...
                repelem(irep, ispecies)'...
                max_xcorrs_allele];

            export_bursts = [export_bursts; ispecies isubnet iparam irep length(bursts) mean(bursts_duration) mean(bursts_total_nodes) mean(bursts_concurrent_nodes)];
            export_per_gene_above = [export_per_gene_above; ispecies isubnet iparam irep per_gene_above_threshold'];
            export_OR = [export_OR; ispecies isubnet iparam irep all_odds_ratio'];
        end
        writematrix(export_per_gene_above, strcat(save_directory, sprintf('per_gene_above-species%d.csv', ispecies)));
        writematrix(all_xcorrs_bursts_export, strcat(save_directory, sprintf('all_xcorrs_burst_species%d-rep%d.csv', ispecies, irep)));
        writematrix(max_xcorrs_bursts_export, strcat(save_directory, sprintf('max_xcorrs_burst_species%d-rep%d.csv', ispecies, irep)));
        writematrix(export_OR, strcat(save_directory, sprintf('long_OR-species%d.csv', ispecies)));
        writematrix(export_bursts, strcat(save_directory, sprintf('long_bursts-species%d.csv', ispecies)));
        writematrix(all_xcorrs_allele_export, strcat(save_directory, sprintf('all_xcorrs_allele_species%d-rep%d.csv', ispecies, irep)));
        writematrix(max_xcorrs_allele_export, strcat(save_directory, sprintf('max_xcorrs_allele_species%d-rep%d.csv', ispecies, irep)));

    end
    
end
