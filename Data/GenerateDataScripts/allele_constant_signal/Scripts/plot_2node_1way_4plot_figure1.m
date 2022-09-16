% Analyze the 2 node, 1 way custom network, A -> B
% The traces of gene expression go to Figure 1D

addpath(genpath(pwd))
addpath(genpath('./../Functions'))
addpath(genpath('./../Data'))
addpath(genpath('~/Documents/MATLAB'))

clear
clc

save_path = './../../../../Paper/plots/plots_for_figures/figure1/20210530_2node_1way_possibilites/';
save_data_path = './../../../../Paper/extractedData/20210528_for_paper/2node_1way/mat_files/';

runID = '20210527_vary_off_add_k_for_paper_1way_2node-rep1';
load(sprintf('S_outpar_%s_2_1_1', runID));
save(strcat(save_data_path, sprintf('S_outpar_%s_2_1_1', runID)), 'S_outpar');

c_11 = distinguishable_colors(11);
c_10 = distinguishable_colors(10);

c_4 = [c_11(1,:); 0.627451, 0.05490196, 0.05490196; c_11(2,:); 0.192157, 0.564709, 0.968627];

%Uncorrelated
%Low noise trans
%211
%OR = 1.2
%0.2 burst transmission
x = 3546;
y = x+120;
f=figure;
f.Position = [100 100 800 300];
tiledlayout(4,1)
for k = 211
    p1=nexttile
    stairs(S_outpar{k}(1, x:y), 'Color', c_11(2,:));
    p2=nexttile
    stairs(S_outpar{k}(3, x:y), 'Color', c_11(2,:));
    p3=nexttile
    stairs(S_outpar{k}(2, x:y), 'Color', c_11(1,:));
    p4=nexttile
    stairs(S_outpar{k}(4, x:y), 'Color', c_11(1,:));
    linkaxes([p1 p2 p3 p4],'xy')

end

saveas(gcf, strcat(save_path,'uncorrelated_low_transmission_4plot_param211.svg'), 'svg');

%medium correlation
%medium noise trans?
%214
%OR = 2.68
%0.6 burst transmission
x = 180;
y = x+120;
f=figure;
f.Position = [100 100 800 300];
tiledlayout(4,1)
for k = 214
    p1=nexttile
    stairs(S_outpar{k}(1, x:y), 'Color', c_11(2,:));
    p2=nexttile
    stairs(S_outpar{k}(3, x:y), 'Color', c_11(2,:));
    p3=nexttile
    stairs(S_outpar{k}(2, x:y), 'Color', c_11(1,:));
    p4=nexttile
    stairs(S_outpar{k}(4, x:y), 'Color', c_11(1,:));
    linkaxes([p1 p2 p3 p4],'xy')

end
saveas(gcf, strcat(save_path,'poised_4plot_param214.svg'), 'svg');

%Correlated
%High noise trans
%218
%OR = 22
%0.85 burst transmission
x = 6000;
y = x+120;
f=figure;
f.Position = [100 100 800 300];
tiledlayout(4,1)
for k = 218
    p1=nexttile
    stairs(S_outpar{k}(1, x:y), 'Color', c_11(2,:));
    p2=nexttile
    stairs(S_outpar{k}(3, x:y), 'Color', c_11(2,:));
    p3=nexttile
    stairs(S_outpar{k}(2, x:y), 'Color', c_11(1,:));
    p4=nexttile
    stairs(S_outpar{k}(4, x:y), 'Color', c_11(1,:));
    linkaxes([p1 p2 p3 p4],'xy')

end
saveas(gcf, strcat(save_path,'corr_high_trans_4plot_param218.svg'), 'svg');
