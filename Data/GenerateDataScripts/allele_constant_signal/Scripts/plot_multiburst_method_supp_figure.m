% Plot an example trace for multi-burst calling
% Supplemental figure

addpath(genpath(pwd))
addpath(genpath('./../Functions'))
addpath(genpath('./../Data'))
addpath(genpath('~/Documents/MATLAB'))

clear
clc
save_path = './../../../../Paper/plots/plots_for_figures/supplemental_figures/';

irep = 1;
runID = sprintf('20210524_vary_off_add_k_for_paper-rep%d', irep);
ispecies = 5;
isubnet= 1;

load(sprintf('S_outpar_%s_%d_%d_1', runID, ispecies, isubnet));
c_11 = distinguishable_colors(11);

summed = cell(1,700);
t=[];
for i = 1:700
    for j = 1:5
        t(j,:) = S_outpar{i}(j,:) + S_outpar{i}(j+5,:);
    end
    summed{i} = t;
end


% short burst
t1=2612;
t2=2661;
for k = [303]
    figure
    for j=1:5
        hold on
        
        stairs(summed{k}(j,t1:t2), 'Color', c_11(j,:))
    end
    legend();
            axis([0 50 0 250]);

end

threshold = 3;
n_alleles = 2;
n_species = 5;

allele_sum = summed{k}(:,t1:t2);

if length(threshold) ~= n_species*n_alleles
    threshold = repmat(threshold,1,n_species*n_alleles);
end

for j = 1:n_species
    above_thresh(j,:) = allele_sum(j,:) > threshold(j);
end
num_nodes_above = sum(above_thresh);

figure
stairs(num_nodes_above, 'Color', 'black');
saveas(gcf, strcat(save_path,'num_nodes_above_short_burst.svg'), 'svg');