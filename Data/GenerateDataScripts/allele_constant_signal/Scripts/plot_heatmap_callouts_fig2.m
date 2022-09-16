%Plot examples for the heatmap in Figure 2B
%Will use rep1, network 5_1

addpath(genpath(pwd))
addpath(genpath('./../Functions'))
addpath(genpath('./../Data'))
addpath(genpath('~/Documents/MATLAB'))

save_path = './../../../../Paper/plots/plots_for_figures/figure2/';
save_data_path = './../../../../Paper/extractedData/20210528_for_paper/vary_off_add_k_700params/mat_files/';

irep = 1;
runID = sprintf('20210604_vary_off_add_k_for_paper_long-rep%d', irep);

ispecies = 5;
isubnet = 1;

load(sprintf('S_outpar_%s_%d_%d_1', runID, ispecies, isubnet));

params_to_plot = [332 334 335 336];

export_data = cell(1,700);

for i = 1:700
    if ismember(i, params_to_plot)
        export_data{i} = S_outpar{i};
    else
        export_data{i} = [0];
    end
end
save(strcat(save_data_path, sprintf('S_outpar_%s_%d_%d_1', runID, ispecies, isubnet)), 'export_data');

c_11 = distinguishable_colors(11);
c_10 = distinguishable_colors(10);

x = 6000;
y = x+500;
for k = 332
    figure
ha(1) = subplot(2,1,1);
for i = 1:5
    hold on
    stairs(S_outpar{k}(i, x:y), 'Color', c_11(i,:));
end
axis([0 500 0 100]);
ha(2) = subplot(2,1,2);
j = 1
for i = 6:10
    hold on
    stairs(S_outpar{k}(i, x:y), 'Color', c_11(j,:));
        j = j+1;
end
axis([0 500 0 100]);

linkaxes(ha,'xy');
end
saveas(gcf, strcat(save_path,'example_OR1_param332.svg'), 'svg');


current = S_outpar{k}(:,100:end);
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

ctable_1 = crosstab(individual_allele_above(1,:), individual_allele_above(6,:));
ctable_4 = crosstab(individual_allele_above(4,:), individual_allele_above(9,:));

save(strcat(save_path,'example_OR1_param332_ctable_node1.mat'), 'ctable_1');
save(strcat(save_path,'example_OR1_param332_all_OR.mat'), 'OR_temp');

x = 6000;
y = x+500;
for k = 334
    figure
ha(1) = subplot(2,1,1);
for i = 1:5
    hold on
    stairs(S_outpar{k}(i, x:y), 'Color', c_11(i,:));
end
axis([0 500 0 100]);

ha(2) = subplot(2,1,2);
j = 1
for i = 6:10
    hold on
    stairs(S_outpar{k}(i, x:y), 'Color', c_11(j,:));
        j = j+1;
end
axis([0 500 0 100]);

linkaxes(ha,'xy');
legend();
end
saveas(gcf, strcat(save_path,'example_OR2_param334.svg'), 'svg');

current = S_outpar{k}(:,100:end);
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
ctable_1 = crosstab(individual_allele_above(1,:), individual_allele_above(6,:));
save(strcat(save_path,'example_OR2_param334_ctable_node1.mat'), 'ctable_1');
save(strcat(save_path,'example_OR2_param334_all_OR.mat'), 'OR_temp');


x = 11000;
y = x+500;
for k = 335
    figure
ha(1) = subplot(2,1,1);
for i = 1:5
    hold on
    stairs(S_outpar{k}(i, x:y), 'Color', c_11(i,:));
end
axis([0 500 0 100]);

ha(2) = subplot(2,1,2);
j = 1
for i = 6:10
    hold on
    stairs(S_outpar{k}(i, x:y), 'Color', c_11(j,:));
        j = j+1;
end
axis([0 500 0 100]);

linkaxes(ha,'xy');
legend();
end
saveas(gcf, strcat(save_path,'example_OR3_param335.svg'), 'svg');

current = S_outpar{k}(:,100:end);
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
ctable_1 = crosstab(individual_allele_above(1,:), individual_allele_above(6,:));
save(strcat(save_path,'example_OR3_param335_ctable_node1.mat'), 'ctable_1');
save(strcat(save_path,'example_OR3_param335_all_OR.mat'), 'OR_temp');

x = 6000;
y = x+500;
for k = 336
    figure
ha(1) = subplot(2,1,1);
for i = 1:5
    hold on
    stairs(S_outpar{k}(i, x:y), 'Color', c_11(i,:));
end
axis([0 500 0 100]);

ha(2) = subplot(2,1,2);
j = 1
for i = 6:10
    hold on
    stairs(S_outpar{k}(i, x:y), 'Color', c_11(j,:));
        j = j+1;
end
axis([0 500 0 100]);

linkaxes(ha,'xy');
legend();
end
saveas(gcf, strcat(save_path,'example_OR5_param336.svg'), 'svg');

current = S_outpar{k}(:,100:end);
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
ctable_1 = crosstab(individual_allele_above(1,:), individual_allele_above(6,:));
save(strcat(save_path,'example_OR5_param336_ctable_node1.mat'), 'ctable_1');
save(strcat(save_path,'example_OR5_param336_all_OR.mat'), 'OR_temp');