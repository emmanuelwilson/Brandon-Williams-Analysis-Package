%%Finds the peak of each cell footprint and compares their proximity
%INPUT: 
% - ms: takse the Miniscope structure
%OUTPUT
% - ndist = returns an n by n matrix of distance between cells where n is the cell number. 
function [ ndist ] = NearestNbr(ms)
ndist = NaN(length(ms.SFPs(1,1,:)),length(ms.SFPs(1,1,:)));

for i = 1: length(ms.SFPs(1,1,:))
    max1 = (max(max(ms.SFPs(:,:,i))));
    [y1, x1] = find((ms.SFPs(:,:,i))== max1);
    for j = 1: length(ms.SFPs(1,1,:))
        max2 = (max(max(ms.SFPs(:,:,j))));
        [y2, x2] = find((ms.SFPs(:,:,j))== max2);
        if j == i
            break
        else 
            ndist(i,j) = sqrt((x2-x1)^2+(y2-y1)^2);                    
        end
    end
end

end