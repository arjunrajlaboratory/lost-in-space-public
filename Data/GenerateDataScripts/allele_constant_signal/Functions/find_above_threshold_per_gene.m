% This function calculates the amount of simulation time that each individual gene is above a user defined expression threshold
% This is used in the analysis of linear networks only

% INPUT:
%
% n_species:    number of species
% n_alleles:    number of alleles (only tested for 2)
% threshold:    threshold for expression
% trim:         trim length to account for initialization time
% S_outpar:     simulation dataset input
%
% OUTPUT:
%
% per_gene_above_threshold:     vector of time above the threshold for each
%                               gene
%
%%
function [per_gene_above_threshold] = ...
    find_above_threshold_per_gene(n_species,n_alleles, threshold,S_outpar, trim)

S_outpar=S_outpar(:,trim:end);

if length(threshold) ~= n_species*n_alleles
    threshold = repmat(threshold,1,n_species*n_alleles);
end

for i = 1:n_species
    allele_sum(i, :) = S_outpar(i,:) + S_outpar(i+n_species,:);
end

for j = 1:n_species
    above_thresh(j,:) = allele_sum(j,:) > threshold(j);
end

per_gene_above_threshold=sum(above_thresh,2);
end
