%Plot examples of increasing input-output correlation
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

%Assessing different params
t1=1;
t2=9000;
for k = [301 302 303 304]
    figure
    for j=1:5
        hold on
        
        stairs(summed{k}(j,t1:t2), 'Color', c_11(j,:))
    end
end
%axis([0 160 0 1]);


%param 301

%single green burst
t1=1210;
t2=1231;
for k = [301]
    figure
    for j=1:5
        hold on
        
        stairs(summed{k}(j,t1:t2), 'Color', c_11(j,:))
    end
                axis([0 25 0 100]);

end
saveas(gcf, strcat(save_path,'train_example_param301_1.svg'), 'svg');

%double, green-red
t1=4131;
t2=4152;
for k = [301]
    figure
    for j=1:5
        hold on
        
        stairs(summed{k}(j,t1:t2), 'Color', c_11(j,:))
    end
                axis([0 25 0 100]);

end
saveas(gcf, strcat(save_path,'train_example_param301_2.svg'), 'svg');

%three, black-green-red
t1=1152;
t2=1173;
for k = [301]
    figure
    for j=1:5
        hold on
        
        stairs(summed{k}(j,t1:t2), 'Color', c_11(j,:))
    end
    legend();
                axis([0 25 0 100]);

end
saveas(gcf, strcat(save_path,'train_example_param301_3.svg'), 'svg');

%five
t1=4451;
t2=4472;
for k = [301]
    figure
    for j=1:5
        hold on
        
        stairs(summed{k}(j,t1:t2), 'Color', c_11(j,:))
    end
    legend();
            axis([0 25 0 100]);

end
saveas(gcf, strcat(save_path,'train_example_param301_4.svg'), 'svg');

% Param 302
%3 burst
t1=3780;
t2=3801;
for k = [302]
    figure
    for j=1:5
        hold on
        
        stairs(summed{k}(j,t1:t2), 'Color', c_11(j,:))
    end
    legend();
            axis([0 25 0 120]);

end
saveas(gcf, strcat(save_path,'train_example_param302_1.svg'), 'svg');

%2 + 3 burst
t1=4264;
t2=4475;
for k = [302]
    figure
    for j=1:5
        hold on
        
        stairs(summed{k}(j,t1:t2), 'Color', c_11(j,:))
    end
    legend();
            axis([0 25 0 120]);

end
saveas(gcf, strcat(save_path,'train_example_param302_2.svg'), 'svg');

%4 burst
t1=128;
t2=149;
for k = [302]
    figure
    for j=1:5
        hold on
        
        stairs(summed{k}(j,t1:t2), 'Color', c_11(j,:))
    end
    legend();
            axis([0 25 0 120]);

end
saveas(gcf, strcat(save_path,'train_example_param302_3.svg'), 'svg');

%4 long
t1=837;
t2=858;
for k = [302]
    figure
    for j=1:5
        hold on
        
        stairs(summed{k}(j,t1:t2), 'Color', c_11(j,:))
    end
    legend();
            axis([0 25 0 120]);

end
saveas(gcf, strcat(save_path,'train_example_param302_4.svg'), 'svg');

% Param 303
% two small burst
t1=380;
t2=429;
for k = [303]
    figure
    for j=1:5
        hold on
        
        stairs(summed{k}(j,t1:t2), 'Color', c_11(j,:))
    end
    legend();
            axis([0 50 0 250]);

end
saveas(gcf, strcat(save_path,'train_example_param303_1.svg'), 'svg');

%two short bursts
t1=825;
t2=874;
for k = [303]
    figure
    for j=1:5
        hold on
        
        stairs(summed{k}(j,t1:t2), 'Color', c_11(j,:))
    end
    legend();
            axis([0 50 0 250]);

end
saveas(gcf, strcat(save_path,'train_example_param303_2.svg'), 'svg');

%two medium bursts
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
saveas(gcf, strcat(save_path,'train_example_param303_3.svg'), 'svg');

% medium burst
t1=5550;
t2=5599;
for k = [303]
    figure
    for j=1:5
        hold on
        
        stairs(summed{k}(j,t1:t2), 'Color', c_11(j,:))
    end
    legend();
            axis([0 50 0 250]);

end
saveas(gcf, strcat(save_path,'train_example_param303_4.svg'), 'svg');

% Param 304

t1=2007;
t2=t1+201;
for k = [304]
    figure
    for j=1:5
        hold on
        
        stairs(summed{k}(j,t1:t2), 'Color', c_11(j,:))
    end
    legend();
        axis([0 200 0 250]);

end
saveas(gcf, strcat(save_path,'train_example_param304_1.svg'), 'svg');

%medium
t1=1878;
t2=t1+201;
for k = [304]
    figure
    for j=1:5
        hold on
        
        stairs(summed{k}(j,t1:t2), 'Color', c_11(j,:))
    end
    legend();
        axis([0 200 0 250]);

end
saveas(gcf, strcat(save_path,'train_example_param304_2.svg'), 'svg');

%smaller, start stop
t1=688;
t2=t1+201;
for k = [304]
    figure
    for j=1:5
        hold on
        
        stairs(summed{k}(j,t1:t2), 'Color', c_11(j,:))
    end
    legend();
        axis([0 200 0 250]);

end
saveas(gcf, strcat(save_path,'train_example_param304_3.svg'), 'svg');

t1=2657;
t2=2858;
for k = [304]
    figure
    for j=1:5
        hold on
        
        stairs(summed{k}(j,t1:t2), 'Color', c_11(j,:))
    end
    legend();
        axis([0 200 0 250]);

end
saveas(gcf, strcat(save_path,'train_example_param304_4.svg'), 'svg');
