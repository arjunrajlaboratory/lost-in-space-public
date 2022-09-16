% Plot example for schematic of signal addition
% Figure 4
addpath(genpath(pwd))
addpath(genpath('./../Functions'))
addpath(genpath('./../Data'))
addpath(genpath('~/Documents/MATLAB'))

save_path = './../../../../Paper/plots/plots_for_figures/figure6/';


load(sprintf('S_outpar_20210518_vary_off_add_k_for_constant_signal-prerep%d_prerun_5_1_1', 1))
prerun_example = cellfun(@(x) [x(1:5,:) + x(6:10,:)], S_outpar, 'UniformOutput', false);
load(sprintf('S_outpar_20210518_vary_off_add_k_for_constant_signal_400-prerep%d_multirep%d_5_1_1', 1, 1))
postrun_exampe = cellfun(@(x) [x(1:5,:) + x(6:10,:); x(11,:)], S_outpar, 'UniformOutput', false);

summed = cell(1,700);
i=1;
for ipre = 1
    for irep = 1:100
        load(sprintf('S_outpar_20210518_vary_off_add_k_for_constant_signal_400-prerep%d_multirep%d_5_1_1', ipre, irep))
        if irep == 1
            summed(i,:) = S_outpar;
        else
            for n = 1:size(summed,2)
                summed{i,n} = summed{i,n} + S_outpar{n};
            end
        end
    end
    i = i+1;
end
mean_summed= cellfun(@(x) x./100, summed, 'UniformOutput', false);
mean_summed_collapse_alleles = cellfun(@(x) x(1:5,:) + x(6:10,:), mean_summed, 'UniformOutput', false);


c_11 = distinguishable_colors(11);

f=figure;
f.Position = [100 100 800 300];
tiledlayout(1,3)
for k= [304]
    nexttile
    for ispecies = 1:5
        hold on
        stairs(prerun_example{1,k}(ispecies,900:end), 'Color', c_11(ispecies,:));
    end
        axis([0 120 0 250]);
legend()
    nexttile
    for ispecies = 1:6
        hold on
        stairs(postrun_exampe{1,k}(ispecies,1:100), 'Color', c_11(ispecies,:));
    end
    %legend();
    axis([0 120 0 259]);
    nexttile
    for ispecies = 1:5
        hold on
                stairs(mean_summed_collapse_alleles{1,k}(ispecies,1:25), 'Color', c_11(ispecies,:));
    end
    axis([0 25 0 250]);

end
saveas(gcf, strcat(save_path,'pre_post_param304.svg'), 'svg');
