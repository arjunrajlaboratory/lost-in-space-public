%Find the cross-correlation between both alleles for each gene in a simulation
%Only tested to work with two n_alleles

%INPUT:
%
% n_species:    number of species
% n_alleles:    number of alleles (only tested for 2)
% trim:         trim length to account for initialization time
% S_outpar:     simulation dataset input
%
%OUTPUT:
%
%all_xcors_export:  cross-correlation for lags up to 10 time units for all n_species
%max_xcorrs_export: maximum cross-correlation and corresponding lag for all n_species
%%
function [all_xcorrs_export, max_xcorrs_export] = ...
    find_xcorrs_alleles(n_species,n_alleles,S_outpar, trim)

S_outpar=S_outpar(:,trim:end);
all_xcorrs_export=[];
max_xcorrs_export=[];

for i = 1:n_species
        [c, lags] = xcorr(S_outpar(i,:), S_outpar(i+n_species,:), 10, 'coeff');
        all_xcorrs_export = [all_xcorrs_export; i c lags];
        [value, index] = max(c);
        max_xcorrs_export = [max_xcorrs_export; i value lags(index)];
end
end