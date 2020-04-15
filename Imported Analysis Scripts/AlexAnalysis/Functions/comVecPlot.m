function params = comVecPlot(data,root)
    
    figure
    set(gcf,'position',[50 50 200.*length(data(:,1)) 200.*length(data(1,:))])

    for i = 1:length(data(:,1))
        for j = i+1:length(data(:,1))
            subplot(length(data(:,1)),length(data(1,:)),(i-1).*length(data(1,:))+j);
            plot(data{i,j}(:,[1 3])',data{i,j}(:,[2 4])','linestyle','-', ...
                'color',[0.75 0.75 0.75],'linewidth',0.5);
            hold on
            plot(data{i,j}(:,[3])',data{i,j}(:,[4])','linestyle','none','color',[0.1 0.9 0.9],...
                'marker','o','markersize',2,'markerfacecolor',[0.1 0.9 0.9])
            set(gca,'xlim',[0 16],'ylim',[0 16])
            axis square
%             axis off
        end
    end
    saveFig(gcf,root,[{'pdf'} {'tiff'}]);
end