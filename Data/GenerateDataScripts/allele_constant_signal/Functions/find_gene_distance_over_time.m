% This function calculates a series of distances for an input simulation
% INPUT:
%
% n_species:    number of species
% n_alleles:    number of alleles (only tested for 2)
% threshold:    threshold for expression
% trim:         trim length to account for initialization time
% S_outpar:     simulation dataset input
%
%
% OUTPUT:
%
% squared_L2_norms_from_median:     matrix of squared L2 norms from median
%                                   expression value of each gene
% squared_L2_norms_from_mean:        matrix of squared L2 norms from mean
%                                   expression value of each gene
% squared_L2_norms_from_zero:       matrix of squared values of gene
%                                   expression
%
%%
function [squared_L2_norms_from_median, squared_L2_norms_from_mean, squared_L2_norms_from_zero] = ...
     find_gene_distance_over_time(n_species,n_alleles,S_outpar, trim)

for i = 1:n_species
    allele_sum(i, :) = S_outpar(i,:) + S_outpar(i+n_species,:);
end

gene_means = mean(allele_sum(:,trim:end),2);
gene_medians = median(allele_sum(:,trim:end),2);

for i = 1:n_species
    squared_L2_norms_from_median(i,:) = (allele_sum(i,:) - gene_medians(i)).^2;
end

squared_L2_norms_from_median(i+1,:) = sum(squared_L2_norms_from_median,1);

for i = 1:n_species
    squared_L2_norms_from_mean(i,:) = (allele_sum(i,:) - gene_means(i)).^2;
end

squared_L2_norms_from_mean(i+1,:) = sum(squared_L2_norms_from_mean,1);


for i = 1:n_species
    squared_L2_norms_from_zero(i,:) = (allele_sum(i,:)).^2;
end

squared_L2_norms_from_zero(i+1,:) = sum(squared_L2_norms_from_zero,1);


end