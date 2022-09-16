%Plot examples of how the OR is calculated for Figure 2A

addpath(genpath(pwd))
addpath(genpath('./../Functions'))
addpath(genpath('./../Data'))
addpath(genpath('~/Documents/MATLAB'))

save_path = './../../../../Paper/plots/plots_for_figures/figure2/';
save_data_path = './../../../../Paper/extractedData/20210528_for_paper/vary_off_add_k_700params/mat_files/';

irep = 1;
runID = sprintf('20210524_vary_off_add_k_for_paper-rep%d', irep);

ispecies = 5;
isubnet = 1;

low_OR_param = 301;

high_OR_param = 304;

load(sprintf('S_outpar_%s_%d_%d_1', runID, ispecies, isubnet));
save(strcat(save_data_path, sprintf('S_outpar_%s_%d_%d_1', runID, ispecies, isubnet)), 'S_outpar');

c_11 = distinguishable_colors(11);

x = 1900;
y = x+500;
for k = 301
    figure
ha(1) = subplot(2,1,1);
for i = 2
    hold on
    stairs(S_outpar{k}(i, x:y), 'Color', c_11(i,:));
end
axis([0 500 0 150]);

ha(2) = subplot(2,1,2);
j = 2
for i = 7
    hold on
    stairs(S_outpar{k}(i, x:y), 'Color', c_11(j,:));
        j = j+1;
end

linkaxes(ha,'x');
axis([0 500 0 150]);
legend();
end
saveas(gcf, strcat(save_path,'example_for_OR_calculation_lowOR_param301.svg'), 'svg');


current = S_outpar{low_OR_param}(:,x:y);
threshold = 3;

threshold = repmat(threshold,1,10);

   individual_allele_above=[];
for i = 1:10
    individual_allele_above(i, :) = current(i,:) > threshold(i);
end



OR_temp = nan(5,1);
for j = 1:5
    temp_table = crosstab(individual_allele_above(j,:), individual_allele_above(j+5,:));
                   if size(temp_table) == [2,2]
                        [h,p,stats] = fishertest(temp_table);
                        OR_temp(j) = stats.OddsRatio;
                   else
                       OR_temp(j) = nan;
                   end
end

node_2_ctable = crosstab(individual_allele_above(2,:), individual_allele_above(7,:));
all_odds_ratio = OR_temp;

save(strcat(save_path,'example_for_OR_calculation_param301_all_OR.mat'), 'all_odds_ratio');
save(strcat(save_path,'example_for_OR_calculation_param301_ctable.mat'), 'node_2_ctable');


x = 1900;
y = x+500;
for k = high_OR_param
    figure
ha(1) = subplot(2,1,1);
for i = 2
    hold on
    stairs(S_outpar{k}(i, x:y), 'Color', c_11(i,:));
end
axis([0 500 0 150]);

ha(2) = subplot(2,1,2);
j = 2
for i = 7
    hold on
    stairs(S_outpar{k}(i, x:y), 'Color', c_11(j,:));
        j = j+1;
end

linkaxes(ha,'x');
axis([0 500 0 150]);
legend();
end
saveas(gcf, strcat(save_path,'example_for_OR_calculation_highOR_param304.svg'), 'svg');

current = S_outpar{high_OR_param}(:,x:y);
threshold = 3;

threshold = repmat(threshold,1,10);

   individual_allele_above=[];
for i = 1:10
    individual_allele_above(i, :) = current(i,:) > threshold(i);
end

OR_temp = nan(5,1);
for j = 1:5
    temp_table = crosstab(individual_allele_above(j,:), individual_allele_above(j+5,:));
                   if size(temp_table) == [2,2]
                        [h,p,stats] = fishertest(temp_table);
                        OR_temp(j) = stats.OddsRatio;
                   else
                       OR_temp(j) = nan;
                   end
end

node_2_ctable = crosstab(individual_allele_above(2,:), individual_allele_above(7,:));
all_odds_ratio = OR_temp;
save(strcat(save_path,'example_for_OR_calculation_param304_all_OR.mat'), 'all_odds_ratio');
save(strcat(save_path,'example_for_OR_calculation_param304_ctable.mat'), 'node_2_ctable');
