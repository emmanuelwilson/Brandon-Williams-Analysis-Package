%Sorts and Rasterizes calcium traces from ms 
function [raster] = raster(ms)
raster = ms.frame;

[~,maxind] = max(raster,[],2);
[~, rastsort] = sort(maxind);
raster = raster(rastsort,:);
end