%Separately generate simplest cases of linear networks of size 2-8 nodes

addpath(genpath(pwd))
addpath(genpath('./../Functions'))
addpath(genpath('./../Data'))

clearvars;
clc;
for ispecies = 2:8
    clearvars -EXCEPT ispecies
    n_species = ispecies;
    
    M = zeros(1,n_species);
    
    for icount = 1:n_species
        ind= [];
        for iM = 1:size(M,1)
            if sum(M(iM,:)) == icount - 1
                ind = [ind,iM];
            end
        end
        for jind = ind
            for icol = 1:n_species
                M_in = M(jind,:);
                M_in(icol) = 1;
                M(length(M)+1,:) = M_in;
            end
        end
        M = unique(M,'rows');
    end
    
    for j = 1:size(M,1)
        M_full{j}(1,:) = M(j,:);
        for ishift = 1:n_species -1
            M_full{j}(ishift+1,:) = circshift(M(j,:),ishift);
        end
    end
    
    count = 1;
    for ifull = 1:length(M_full)
        G = digraph(M_full{ifull});
        weak_bins = conncomp(G,'Type','weak');
        if all(weak_bins == 1)
            M_full_con{count} = M_full{ifull};
            count = count + 1;
        end
    end
    
    
    for ifullunique = 1:length(M_full_con)
        if ifullunique == 1
            M_iso{1} = M_full_con{1};
        else
            G1 = digraph(M_full_con{ifullunique});
            countiso = 0;
            for iiso = 1:length(M_iso)
                G2 = digraph(M_iso{iiso});
                if isisomorphic(G1, G2) == 1
                    break;
                else
                    countiso = countiso + 1;
                end
            end
            if countiso == length(M_iso)
                M_iso{length(M_iso)+1} = M_full_con{ifullunique};
            end
        end
    end
    
    M_iso = M_iso{1};
    M_iso(1,:) = repelem(0, ispecies);
    M_iso = {M_iso};
    %M_iso_save = sprintf('./../Data/M_iso_linear%d', ispecies);
    %save(M_iso_save,'M_iso')
end