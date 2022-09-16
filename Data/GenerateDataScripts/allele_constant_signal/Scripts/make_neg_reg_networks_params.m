%Generate unique parameter ensemble and network for negative regulation
%tests

addpath(genpath(pwd))
addpath(genpath('./../Functions'))
addpath(genpath('./../Data'))
clear
clc
%Generate parameter set with large basal on rates to allow for subtraction
%of r_add

r_deg1_to_test = .01;
r_prod1_ratios = .01;
r_off1_ratios = logspace(0, 2, 10);
r_add1_ratios = .025;
r_on1_ratios = logspace(-1, 2, 10);

proddiff1 = 10000;
k1 = linspace(20, 200, 7);

n1 = 1;
p = 1;
combos = combvec(r_prod1_ratios,...
    r_deg1_to_test, r_add1_ratios, ...
    r_off1_ratios, proddiff1,...
    r_on1_ratios, k1, n1, p)';

corrected_combos = [combos(:,1).*combos(:,2) combos(:,2) combos(:,3).*combos(:,2).*combos(:,9)...
    combos(:,4).*combos(:,2).*combos(:,9) combos(:,5) combos(:,6).*combos(:,2).*combos(:,9) combos(:,7:8) combos(:,9)];
DataParams = corrected_combos;
save('./../Data/20220318_negative_regulation_params','DataParams')

%Generate example 5-node representative networks that include negative
%regulation

conn_1_all_negative = {[0 0 0 0 -1; -1 0 0 0 0; 0 -1 0 0 0; 0 0 -1 0 0; 0 0 0 -1 0]};
conn_4_one_negative = {[0 1 1 1 -1; -1 0 1 1 1; 1 -1 0 1 1; 1 1 -1 0 1; 1 1 1 -1 0]};
conn_4_all_negative = {[0 -1 -1 -1 -1; -1 0 -1 -1 -1; -1 -1 0 -1 -1; -1 -1 -1 0 -1; -1 -1 -1 -1 0]};

M_iso = [conn_1_all_negative, conn_4_one_negative, conn_4_all_negative];

M_iso_save = sprintf('./../Data/M_iso_neg_reg%d', 5);
save(M_iso_save,'M_iso')