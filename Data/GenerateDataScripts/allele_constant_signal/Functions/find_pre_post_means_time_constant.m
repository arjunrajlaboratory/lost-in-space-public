function[mean_value_of_pre_signal_collapse_alleles, mean_value_of_post_signal_replicates_collapse_alleles, time_constant] = ...
    find_pre_post_means_time_constant(n_species,n_alleles, subnet, run_id, num_pre_replicates, num_post_replicates, nparams, trim, prerunID)
summed_over_post_signal_replicates = cell(1,nparams);
i=1;
for ipre = 1:num_pre_replicates
    for irep = 1:num_post_replicates
        load(sprintf('S_outpar_%s-prerep%d_multirep%d_%d_%d_1', run_id, ipre, irep, n_species, subnet))
        if irep == 1
            summed_over_post_signal_replicates(i,:) = S_outpar;
        else
            for n = 1:size(summed_over_post_signal_replicates,2)
                summed_over_post_signal_replicates{i,n} = summed_over_post_signal_replicates{i,n} + S_outpar{n};
            end
        end
    end
    i = i+1;
end
mean_of_post_signal_replicates= cellfun(@(x) x./num_post_replicates, summed_over_post_signal_replicates, 'UniformOutput', false);
if n_alleles == 2
    mean_of_post_signal_replicates_collapse_alleles = cellfun(@(x) x(1:n_species,:) + x((n_species+1):(n_species*2),:), mean_of_post_signal_replicates, 'UniformOutput', false);
else
    error('only works for 2 alleles for now!');
end
mean_value_of_post_signal_replicates_collapse_alleles = cellfun(@(x) mean(x, 2), mean_of_post_signal_replicates_collapse_alleles, 'UniformOutput', false);

time_constant = zeros(n_species,num_pre_replicates,nparams);
for ipre = 1:num_pre_replicates
    for j = 1:nparams
        for i=1:n_species

            yi=y(end);
            idx=max(find(abs(y-yi)>=0.37*yi));

            f = find(diff([0,mean_of_post_signal_replicates_collapse_alleles{ipre,j}(i,2:end) <= mean_value_of_post_signal_replicates_collapse_alleles{ipre,j}(i,:).*0.632,0]==1));
            p = f(1:2:end-1);  
            y = f(2:2:end)-p;  
            idx = find(p==1);
            
            if ~isempty(idx)
                time_constant(i,ipre,j) = y(idx);
            else
                time_constant(i,ipre,j) = NaN;
            end
            
        end
    end
end

mean_value_of_pre_signal_collapse_alleles = cell(num_pre_replicates,nparams);
for ipre = 1:num_pre_replicates
    load(sprintf('S_outpar_%s%d_prerun_%d_%d_1', prerunID, ipre, n_species, subnet));
    tcollapse = cellfun(@(x) x(1:n_species,:) + x((n_species+1):(n_species*2),:), S_outpar, 'UniformOutput', false);
    mean_value_of_pre_signal_collapse_alleles(ipre,:)=cellfun(@(x) mean(x(:,trim:end), 2), tcollapse, 'UniformOutput', false);
end

end