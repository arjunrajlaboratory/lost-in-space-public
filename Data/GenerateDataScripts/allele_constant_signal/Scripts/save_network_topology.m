%Extract characteristic distance, connectivity, self-loops for all networks

addpath(genpath(pwd))
addpath(genpath('./../Functions'))
addpath(genpath('./../Data'))

save_directory = './../../../../Paper/extractedData/20210528_for_paper/';

count = 1;
for n_species = 2:8
    
    load(sprintf('M_iso%d', n_species));
    
    for inet = 1:length(M_iso)
    
        n = size(M_iso{inet},1);
        
        for in = 1:n
            for jn = 1:n
                D{count}{inet}(in,jn) = length(shortestpath(digraph(M_iso{inet}),in,jn))-1;
            end
        end
        
        for in= 1:n
            if M_iso{inet}(in,in) == 1
                D{count}{inet}(in,in) = D{count}{inet}(in,in)+1;
            else
                ind = find(M_iso{inet}(:,in)==1);
                I = [];
                for iind = ind'
                    I = [I,length(shortestpath(digraph(M_iso{inet}),in,iind))-1];
                end
                D{count}{inet}(in,in) = D{count}{inet}(in,in) + min(I) + 1;
            end
        end

        L{count}(inet) = mean(mean(D{count}{inet}))/n;
        conn{count}(inet) = sum(M_iso{inet}(1,:));
        
        loops{count}(inet) =  M_iso{inet}(1,1);
    end
    count = count + 1;
end

count = 1;
topology_export=[];
for n_species = 2:8
   for isubnet = 1:length(loops{count})
       topology_export = [topology_export; n_species, isubnet, conn{count}(isubnet), loops{count}(isubnet), L{count}(isubnet)];
   end
    count = count + 1;
end

writematrix(topology_export, strcat(save_directory, 'topology_metadata.csv'))