
function help_showpairwise(m)
    if ~iscell(m)
        m = num2cell(m);
    end

    figure
    step = 0.01;
    cm = [[[0:step:1]'; ones(length(0:step:1),1)] ...
        [[0:step:1 1:-step:0]'] ...
        flipud([[0:step:1]'; ones(length(0:step:1),1)])];
%         imagesc(crossSim)

    toPlot = cellfun(@nanmedian,m);
    toPlot = nanmax(toPlot,toPlot');
%     order = [1 7 6 3 4 5 2 8];
    order = 1:length(toPlot);
    toPlot = toPlot(order,order);
    imagesc(toPlot)
    colorbar('color','w')
%     caxis([0.2 0.8])
    alpha(double(~isnan(toPlot)))
    xlabel('Session','color','w')
    ylabel('Session','color','w')
%     colormap(cm);
    colormap('parula')
    set(gcf,'color','k')
    set(gca,'color','k','xgrid','off','ygrid','off')
end





