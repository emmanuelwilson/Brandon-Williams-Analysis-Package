%% Will sort your calcium traces by chronological maximal peak
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   INPUT:  
%       -ms: miniscope structure:
%           *ms.FiltTraces, Calcium transient for everycell
%   OUTPUT:
%       -raster: m by n matrix of your reorganized calcium traces. m
%       represents the cell and n time. 
%       -rastsort: Index used to create reorganized calcium traces ie.
%       raster = A(rastsort,:)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [raster, rastsort] = ms_sortpeaks(ms)
ms.FiltTraces = rot90(ms.FiltTraces);
[~,maxind] = max(ms.FiltTraces,[],2);
[~, rastsort] = sort(maxind);
raster = ms.FiltTraces(rastsort,:);
end