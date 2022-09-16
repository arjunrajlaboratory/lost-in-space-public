%A revised version of find_RNA_bursts function that allows a user-defined
%'burst_window" to call multi-bursts, allowing a short peroid of time in
%which no gene is above the threshold to still be called a multi-burst
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
%filter_burst_bool: boolean whether to require that 'multi-bursts' contain
%               at least two different genes to surpass the threshold to
%               count as a multi-burst. TRUE is used for all analyses in
%               the paper
%burst_window:  the size of the window for a moving average with which to
%               find multi-bursts
%OUTPUT:
%
%bursts:       all multi-bursts in the simulation
%bursts_duration:   vector of how long each multi-burst lasts
%bursts_total_nodes:    vector of the total genes that surpassed the
%                       threshold during each multi-burst
%bursts_concurrent_nodes:   vector of the total genes that concurrently
%                           surpassed the threshold during each multi-burst
%nodes_fraction_on:         matrix containing the fraction of time each
%                           gene was above the threshold during each multi-burst
%
%%
function [bursts, bursts_duration, bursts_total_nodes, bursts_concurrent_nodes, nodes_fraction_on] = ...
     find_RNA_bursts_with_window(n_species,n_alleles, threshold,S_outpar, trim, filter_burst_bool, burst_window)

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

rolling_above = movsum(above_thresh,burst_window,2)>0;
num_nodes_above = sum(rolling_above);

threshold2 = 0;
transitions = diff([0, num_nodes_above > threshold2, 0]);
runstarts = find(transitions == 1);
runends = find(transitions == -1) - 1;
tbursts = arrayfun(@(s, e) rolling_above(:,s:e), runstarts, runends, 'UniformOutput', false);
if filter_burst_bool
    idx = cellfun(@(x) sum(any(x,2)) > 1, tbursts);
    bursts = tbursts(idx);
else
    bursts = tbursts;
end
bursts_duration = zeros(length(bursts), 1);
bursts_total_nodes = zeros(length(bursts), 1);
bursts_concurrent_nodes = zeros(length(bursts), 1);
nodes_fraction_on = zeros(length(bursts), n_species);
for i = 1:length(bursts)
    bursts_duration(i)= size(bursts{i},2);
    bursts_total_nodes(i) = sum((sum(bursts{i},2)>0));
    bursts_concurrent_nodes(i) = max(sum(bursts{i},1));
    nodes_fraction_on(i,:) = sum(bursts{i},2)/size(bursts{i},2);
end

end