%This function measures allele correlation metrics for a given simuation
%output. Binarizes allele level data and calculates odds ratio across the
%simulation for each node
%
%
%INPUT:
%
%n_species:     number of nodes in the network of given simulation (eg. 2, 3, 5, 8)
%n_alleles:     number of alleles in the network of given simulation
%threshold:     threshold above which an allele is said to be expressing
%               for binarization purposes
%S_outpar:      output of simulation file - number of gene product counts per
%               time unit per gene and state of DNA per gene per time unit (on/off)
%trim:          amount to trim from the start of a simulation to account for initialization time 
%OUTPUT:
%
%mean_odds_ratio:       mean of binarized odds ratio of allelic correlation
%                       over all nodes
%all_odds_ratio:        all of the odds ratios, not averaged
%%
function [mean_odds_ratio, all_odds_ratio]...
    = find_RNA_allele_OR(n_species,n_alleles,threshold,S_outpar,trim)

S_outpar = S_outpar(:,trim:end);

if length(threshold) ~= n_species*n_alleles
   threshold = repmat(threshold,1,n_species*n_alleles);
end

if n_alleles > 2
    error("only works for 2 allele now");
end
individual_allele_above=[];
for i = 1:n_species*n_alleles
    individual_allele_above(i, :) = S_outpar(i,:) > threshold(i);
end


OR_temp = nan(n_species,1);
for j = 1:n_species
    temp_table = crosstab(individual_allele_above(j,:), individual_allele_above(j+n_species,:));
                   if size(temp_table) == [2,2]
                        [h,p,stats] = fishertest(temp_table);
                        OR_temp(j) = stats.OddsRatio;
                   else
                       OR_temp(j) = nan;
                   end
end
all_odds_ratio = OR_temp;
mean_odds_ratio = mean(OR_temp,'omitnan');
end