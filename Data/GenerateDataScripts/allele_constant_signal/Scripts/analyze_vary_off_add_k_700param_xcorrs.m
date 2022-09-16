%Analyze all networks for the vary off add k 700 parameter ensemble using
%cross correlation

addpath(genpath(pwd))
addpath(genpath('./../Functions'))
addpath(genpath('./../Data'))

clear
clc
save_directory = './../../../../Paper/extractedData/20210528_for_paper/cross_correlation/';
threshold = 3;
n_params = 700;
trim_length = 100;
n_alleles = 2;


for irep = 1
    irep
    runID = sprintf('20210604_vary_off_add_k_for_paper_long-rep%d', irep);

    min_species = 2;
    max_species = 5;
    all_xcorrs_export=[];
    max_xcorrs_export=[];
    
    for ispecies = min_species:max_species
        ispecies
        load(sprintf('M_iso%d', ispecies))
        max_subnet = size(M_iso,2);
        export_rows = ispecies * (ispecies - 1) / 2;
        all_xcorrs_bursts_export=[];
        max_xcorrs_bursts_export=[];
        all_xcorrs_allele_export=[];
        max_xcorrs_allele_export=[];

        for isubnet = 1:max_subnet

            load(sprintf('S_outpar_%s_%d_%d_1', runID, ispecies, isubnet));


            for iparam = 1:n_params

                [all_xcorrs_bursts, max_xcorrs_bursts]=...
                    find_xcorrs_RNA_bursts(ispecies,n_alleles,S_outpar{iparam},trim_length);
                [all_xcorrs_allele, max_xcorrs_allele]=...
                    find_xcorrs_alleles(ispecies,n_alleles,S_outpar{iparam},trim_length);

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
                
            end
        end
    writematrix(all_xcorrs_bursts_export, strcat(save_directory, sprintf('all_xcorrs_burst_species%d-rep%d.csv', ispecies, irep)));
    writematrix(max_xcorrs_bursts_export, strcat(save_directory, sprintf('max_xcorrs_burst_species%d-rep%d.csv', ispecies, irep)));
    writematrix(all_xcorrs_allele_export, strcat(save_directory, sprintf('all_xcorrs_allele_species%d-rep%d.csv', ispecies, irep)));
    writematrix(max_xcorrs_allele_export, strcat(save_directory, sprintf('max_xcorrs_allele_species%d-rep%d.csv', ispecies, irep)));

    end
end

load('20210407_vary_off_add_k_fewer_params_for_signalburst')
writematrix(DataParams, strcat(save_directory, 'metadata.csv'));