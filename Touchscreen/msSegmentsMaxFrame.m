function [Scalledfoot] = msSegmentsMaxFrame(ms,segments)%(ms, frames, filtering)
%%Splits activity in the specified # of segments and creates maximal projects for each.
%This function will normalize the calcium trace (from 0 to 1) and scale the
%brightness of the total maximal projection to the local segment maxima.
%Creating a series of spatial footprints representative of the maximal
%projects for each segment. 
%   INPUT:
%           -ms: miniscope structure, need these subparameters
%               *SFPs: Spatial Footprints of each cell. n by m by p, where
%               p is the cell# and n and m the spatial FOV.
%           -segments: User defined number of segments to separate activity
%   OUTPUT:
%           -Scalledfoot: n-by-m-by-p-by-s 4 dimensional matrix. Where n
%           and m are the spatial representation of cell p in segment s.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Scalledfoot = zeros(length(ms.SFPs(:,1,1)),length(ms.SFPs(1,:,1)),length(ms.SFPs(1,1,:)),segments);
vidind = zeros(segments,1);
vidmaxF = zeros(length(ms.FiltTraces(1,:)),segments);
ind = 0;
for i = 1 : segments
    ind = ind + round(length(ms.FiltTraces(:,1))/segments);
    vidind(i) = ind;
end
for i = 1: length(ms.FiltTraces(1,:))
    temp = normalize(ms.FiltTraces(:,i),'range');
    for j = 1 : segments
        if j == 1
        vidmaxF(i,j) = max(temp(1:vidind(j)));
        elseif j == segments
            vidmaxF(i,j) = max(temp(vidind(j-1),end));
        else
            vidmaxF(i,j) = max(temp(vidind(j-1):vidind(j)));
        end
        Scalledfoot(:,:,i,j) = ms.SFPs(:,:,i).*vidmaxF(i,j);
    end
    
end

end

