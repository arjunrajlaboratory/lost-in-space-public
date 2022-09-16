%Generate 3D plots for figure 3

addpath(genpath(pwd))
addpath(genpath('./../Functions'))
addpath(genpath('./../Data'))
addpath(genpath('~/Documents/MATLAB'))

clear
clc
rng(42)


save_path = './../../../../Paper/plots/plots_for_figures/';
save_data_path = './../../../../Paper/extractedData/20210528_for_paper/vary_off_add_k_700params/mat_files/';

az = 100;
elv= 12;

viridis_mako_3 = ["#0D0887" "#CC4678" "#000000"];
irep = 1;
runID = sprintf('20210524_vary_off_add_k_for_paper-rep%d', irep);
ispecies = 3;
isubnet= 1;

load(sprintf('S_outpar_%s_%d_%d_1', runID, ispecies, isubnet));
c_11 = distinguishable_colors(11);

summed = cell(1,700);
t=[];
for i = 1:700
    for j = 1:3
        t(j,:) = S_outpar{i}(j,:) + S_outpar{i}(j+3,:);
    end
    summed{i} = t;
end

save(strcat(save_data_path, sprintf('S_outpar_%s_%d_%d_1_summed.mat', runID, ispecies, isubnet)), 'summed');

n_alleles=2;
trim_length=100;
squared_L2_norms_from_mean=cell(1,700);
squared_L2_norms_from_zero=cell(1,700);
squared_L2_norms_from_median=cell(1,700);
for iparam = 1:700
    [squared_L2_norms_from_median{iparam} squared_L2_norms_from_mean{iparam}, squared_L2_norms_from_zero{iparam}]=...
        find_gene_distance_over_time(ispecies,n_alleles, S_outpar{iparam}, trim_length);
end
save(strcat(save_data_path, sprintf('squared_L2_from_median_%s_%d_%d.mat', runID, ispecies, isubnet)), 'squared_L2_norms_from_median');
save(strcat(save_data_path, sprintf('squared_L2_from_mean_%s_%d_%d.mat', runID, ispecies, isubnet)), 'squared_L2_norms_from_mean');
save(strcat(save_data_path, sprintf('squared_L2_from_zero_%s_%d_%d.mat', runID, ispecies, isubnet)), 'squared_L2_norms_from_zero');

idx = [341 302 304];

t1=1280;
t2=1396;
for k = 1:length(idx)
    figure
    for j=1:3
        hold on
        
        stairs(summed{idx(k)}(j,t1:t2), 'Color', c_11(j,:))
    end
    if k ==1
        ylim([0 35])
    elseif k == 2
        ylim([0 175])
    elseif k == 3
        ylim([0 225])
    end
    saveas(gcf, strcat(save_path,'trace_param', string(k), '.svg'), 'svg');
end
jitterScale = 1.0;
jitterX = rand(size(summed{1},2),1) * jitterScale;
jitterY = rand(size(summed{1},2),1) * jitterScale;
jitterZ = rand(size(summed{1},2),1) * jitterScale;

jitterMatrix = [jitterX, jitterY, jitterZ]';

jitteredData = cellfun(@(x) x(:,:) + jitterMatrix, summed, 'UniformOutput', false);

f=figure;
f.Position = [100 100 800 300];
tiledlayout(1,3)
for k = 1:length(idx)
subplot(1,3,k)
colormap('magma')
patch([jitteredData{idx(k)}(1,t1:t2) nan], [jitteredData{idx(k)}(2,t1:t2) nan], [jitteredData{idx(k)}(3,t1:t2) nan],[1:117 nan], 'FaceColor','none','EdgeColor','interp')
axis equal
hold on
tmedians= median(summed{idx(k)},2);

scatter3(tmedians(1), tmedians(2), tmedians(3), 75, 'fill', 'MarkerFaceColor', 'black');

if k ==1
    xlim([0 35])
    ylim([0 35])
    zlim([0 35])
elseif k == 2 
    xlim([0 175])
    ylim([0 175])
    zlim([0 175])
elseif k == 3
    xlim([0 225])
    ylim([0 225])
    zlim([0 225])
end
grid on
view([az, elv])
% hAxis = gca;
% hAxis.ZAxis.FirstCrossoverValue  = hAxis.XLim(2); 
% hAxis.ZAxis.SecondCrossoverValue = hAxis.XLim(2);
% hAxis.ZAxis.FirstCrossoverValue  = hAxis.XLim(1);
% hAxis.ZAxis.SecondCrossoverValue = hAxis.XLim(1);
% hAxis.XAxis.FirstCrossoverValue  = hAxis.XLim(1);
% hAxis.XAxis.SecondCrossoverValue = hAxis.XLim(1);
% hAxis.YAxis.FirstCrossoverValue  = hAxis.XLim(1);
% hAxis.YAxis.SecondCrossoverValue = hAxis.XLim(1);
end
set(gcf,'Renderer','Painter')
hgexport(gcf,strcat(save_path,'faceted_LIS_params341_304_306.eps'));
colorbar
hgexport(gcf,strcat(save_path,'colorbar.eps'));
