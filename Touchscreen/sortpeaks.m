%Will sort through the raster plot reoganizing cells by time of highest
%peak and return new raster plot
function [rast,rastsort] = sortpeaks(rast)
[~,maxind] = max(rast,[],2);
[~, rastsort] = sort(maxind);
rast = rast(rastsort,:);
end