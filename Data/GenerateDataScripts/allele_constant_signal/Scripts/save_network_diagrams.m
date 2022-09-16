%Save network diagrams as PDF for importing into supplemental figure

clearvars -except saveFigures
clc

load('./../Data/M_iso2')
count = 0;

for j = 1:length(M_iso)
    
    pos = [0.1,0.32,0.54,0.76];
    xlen = 0.2;
    ylen = 0.2;
    
    figure
    i = 1;
    subplot('Position',[pos(i),pos(1),xlen,ylen]);
    
    count = count + 1;
    G = digraph(M_iso{j});
    h = plot(G,'MarkerSize',5, 'NodeLabel',{}, 'LineWidth', 1, 'ArrowSize', 5, 'Edgecolor','k','Nodecolor',[150,150,150]./255,'EdgeAlpha',1);
    ax = gca;
    ax.Visible = 'off';
    

        set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
        path = sprintf('./../Plots/Figure3A%d',j);
        print('-dpdf',[path])

    
end


load('./../Data/M_iso3')
count = 0;

for j = 1:length(M_iso)
    
    pos = [0.1,0.32,0.54,0.76];
    xlen = 0.2;
    ylen = 0.2;
    
    figure
    i = 1;
    subplot('Position',[pos(i),pos(1),xlen,ylen]);
    
    count = count + 1;
    G = digraph(M_iso{j});
    h = plot(G,'MarkerSize',5, 'NodeLabel',{}, 'LineWidth', 1, 'ArrowSize', 5, 'Edgecolor','k','Nodecolor',[150,150,150]./255,'EdgeAlpha',1);
    ax = gca;
    ax.Visible = 'off';
    

        set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
        path = sprintf('./../Plots/Figure3B%d',j);
        print('-dpdf',[path])

    
end

load('./../Data/M_iso4')
count = 0;

for j = 1:length(M_iso)
    
    pos = [0.1,0.32,0.54,0.76];
    xlen = 0.2;
    ylen = 0.2;
    
    figure
    i = 1;
    subplot('Position',[pos(i),pos(1),xlen,ylen]);
    
    count = count + 1;
    G = digraph(M_iso{j});
    h = plot(G,'MarkerSize',5, 'NodeLabel',{}, 'LineWidth', 1, 'ArrowSize', 5, 'Edgecolor','k','Nodecolor',[150,150,150]./255,'EdgeAlpha',1);
    ax = gca;
    ax.Visible = 'off';
    

        set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
        path = sprintf('./../Plots/Figure3C%d',j);
        print('-dpdf',[path])

    
end


load('./../Data/M_iso5')
count = 0;

for j = 1:length(M_iso)
    
    pos = [0.1,0.32,0.54,0.76];
    xlen = 0.2;
    ylen = 0.2;
    
    figure
    i = 1;
    subplot('Position',[pos(i),pos(1),xlen,ylen]);
    
    count = count + 1;
    G = digraph(M_iso{j});
    h = plot(G,'MarkerSize',5, 'NodeLabel',{}, 'LineWidth', 1, 'ArrowSize', 5, 'Edgecolor','k','Nodecolor',[150,150,150]./255,'EdgeAlpha',1);
    ax = gca;
    ax.Visible = 'off';
    

        set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
        path = sprintf('./../Plots/Figure3D%d',j);
        print('-dpdf',[path])

    
end


load('./../Data/M_iso_neg_reg5')
count = 0;

for j = 1:length(M_iso)
    
    pos = [0.1,0.32,0.54,0.76];
    xlen = 0.2;
    ylen = 0.2;
    
    figure
    i = 1;
    subplot('Position',[pos(i),pos(1),xlen,ylen]);
    
    count = count + 1;
    G = digraph(M_iso{j});
    h = plot(G,'MarkerSize',5, 'NodeLabel',{}, 'LineWidth', 1, 'ArrowSize', 5, 'Edgecolor','k','Nodecolor',[150,150,150]./255,'EdgeAlpha',1);
    ax = gca;
    ax.Visible = 'off';
    

        set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
        path = sprintf('./../Plots/Figure3E%d',j);
        print('-dpdf',[path])

    
end



count = 0;
for k = 2:5
    load(sprintf('./../Data/M_iso_linear%d',k))
    
    for j = 1:length(M_iso)
        pos = [0.1,0.32,0.54,0.76];
        xlen = 0.2;
        ylen = 0.2;
        
        figure
        i = 1;
        subplot('Position',[pos(i),pos(1),xlen,ylen]);
        
        count = count + 1;
        G = digraph(M_iso{j});
        h = plot(G,'MarkerSize',5, 'NodeLabel',{}, 'LineWidth', 1, 'ArrowSize', 5, 'Edgecolor','k','Nodecolor',[150,150,150]./255,'EdgeAlpha',1);
        ax = gca;
        ax.Visible = 'off';
        
        
        set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
        path = sprintf('./../Plots/Figure3F%d',k);
        print('-dpdf',[path])
        
        
    end
end
