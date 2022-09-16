% Calculate pre and post signal means and time constant
% Figure 4
addpath(genpath(pwd))
addpath(genpath('./../Functions'))
addpath(genpath('./../Data'))

clear
clc

save_path = './../../../../Paper/extractedData/20210528_for_paper/vary_off_add_k_700params/pre_post_signal/';

%Analayse multireplicate from multiple initialization conditions

%we have initialization replicates 1-100, each with 100 replicates


%node 5, summed over replicates
summed = cell(1,700);
i=1;
for ipre = 1:100
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

starting_values = cellfun(@(x) x(:,1), mean_summed_collapse_alleles, 'UniformOutput', false);
starting_values_export=[];
for ipre = 1:100
    for jparam = 1:700
        starting_values_export = [starting_values_export; ipre jparam starting_values{ipre,jparam}'];
    end
end

writematrix(starting_values_export, strcat(save_path, 'starting_values.csv'));

post_signal_means = cellfun(@(x) mean(x(:,10:end), 2), mean_summed_collapse_alleles, 'UniformOutput', false);

pre_signal_means = cell(100,700);
for ipre = 1:100
    load(sprintf('S_outpar_20210518_vary_off_add_k_for_constant_signal-prerep%d_prerun_5_1_1', ipre));
    tcollapse = cellfun(@(x) x(1:5,:) + x(6:10,:), S_outpar, 'UniformOutput', false);
    pre_signal_means(ipre,:)=cellfun(@(x) mean(x(:,100:end), 2), tcollapse, 'UniformOutput', false);
end

temp=permute(reshape([pre_signal_means{:}], [], size(pre_signal_means, 1), size(pre_signal_means, 2)), [2 3 1]);
grand_pre_signal_means = squeeze(mean(temp,1))';


time_constant = zeros(5,100,700);
for ipre = 1:100
    for j = 1:700
        for i=1:5
            response_value = (post_signal_means{ipre,j}(i,:) - grand_pre_signal_means(i,j)) * 0.666 + grand_pre_signal_means(i,j);
            if (mean_summed_collapse_alleles{ipre,j}(i,1) >= response_value)
                tc = nan;
            else
                f = find(diff([0,mean_summed_collapse_alleles{ipre,j}(i,2:end) >= response_value,0]==1));
                tc = f(1);
            end
            time_constant(i,ipre,j) = tc;
            
        end
    end
end
export_pre_signal_means = [];
for ipre = 1:100
    for jparam = 1:700

            export_pre_signal_means = [export_pre_signal_means; ipre jparam pre_signal_means{ipre,jparam}'];

    end
end

export_time_constant = [];
for ipre = 1:100
    for jparam = 1:700
        for ispecies = 1:5
            export_time_constant = [export_time_constant; ipre jparam ispecies time_constant(ispecies, ipre, jparam)];
        end
    end
end

export_post_signal_means = [];
for ipre = 1:100
    for jparam = 1:700
        export_post_signal_means = [export_post_signal_means; ipre jparam post_signal_means{ipre,jparam}'];
        
    end
end
writematrix(export_pre_signal_means, strcat(save_path, 'pre_signal_means.csv'));
writematrix(export_post_signal_means, strcat(save_path, 'post_signal_means.csv'));
writematrix(export_time_constant, strcat(save_path, 'export_time_constant.csv'));
