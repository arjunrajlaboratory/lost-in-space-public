% This function calculates the fraction of simulation time that all genes are above a user defined expression threshold
% This is used to call whether a given simulation is constitutively high

% INPUT:
%
% n_species:    number of species
% n_alleles:    number of alleles (only tested for 2)
% threshold:    threshold for expression
% trim:         trim length to account for initialization time
% S_outpar:     simulation dataset input
% upper_bound:  the fraction of time above threshold needed to be called
%               constitutively high
%
% OUTPUT:
%
% fraction all _above:     numeric fraction of simulation time that all
%                          genes are above the threshold
% is_constant_high:        boolean for whether simulation is constitutively
%                          expressing according to the upper_bound
%
%%
function [fraction_all_above, is_constant_high] =...
    find_fraction_all_above_threshold(n_species,n_alleles, threshold,S_outpar, trim, upper_bound)

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

all_above = sum(above_thresh, 1) == n_species;

fraction_all_above = sum(all_above)/length(all_above);
is_constant_high = 0;

if fraction_all_above > upper_bound
    is_constant_high = 1;
end
end