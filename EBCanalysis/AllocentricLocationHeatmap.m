%% Will take your tracking and make an allocentric heatmap of mouse location
function [allmap] = AllocentricLocationHeatmap(tracking,context,mouse)
map = AllocentricTraceMaps(tracking(:,1:2)',ones(length(tracking(:,1)))',ones(length(tracking(:,1))));
allmap = map(:,:,1);
allmap(find(isnan(allmap))) = 0;
figure(1)
b = imagesc(allmap)
% set(b,'AlphaData',~isnan(allmap(:,:)))
title(['Context ', context, ' Location Heatmap Mouse: ', mouse],'FontSize',26)
saveas(gcf,['LocationHeatmap_Context',context,'_', mouse, '.fig'])
saveas(gcf,['LocationHeatmap_Context',context,'_', mouse, '.eps'])
end
