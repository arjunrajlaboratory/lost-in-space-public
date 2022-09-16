%For analysis of 2 node 1 way, not ultimately used in paper
function [bursts, bursts_duration] = ...
     find_upstream_bursts(n_species,n_alleles, threshold,S_outpar, trim)

 S_outpar=S_outpar(:,trim:end);
 if length(threshold) ~= n_species*n_alleles
     threshold = repmat(threshold,1,n_species*n_alleles);
 end
 
 for i = 1:n_species
     allele_sum(i, 1:size(S_outpar,2)) = S_outpar(i,:) + S_outpar(i+n_species,:);
 end
 
 for j = 1:n_species
     above_thresh(j,1:size(allele_sum,2)) = allele_sum(j,:) > threshold(j);
 end
 num_nodes_above = sum(above_thresh);
 
 transitions = diff([0, above_thresh(1,:) > 0, 0]);
 runstarts = find(transitions == 1);
 runends = find(transitions == -1) - 1;
 bursts = arrayfun(@(s, e) above_thresh(1,s:e), runstarts, runends, 'UniformOutput', false);
 bursts_duration = zeros(length(bursts), 1);

 for i = 1:length(bursts)
     bursts_duration(i)= size(bursts{i},2);
 end

end