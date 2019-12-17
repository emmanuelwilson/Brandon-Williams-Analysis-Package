%%Compares difference between two values (angle or distance) with pixel
%%distance.
%INPUT:
%- var: variable you want to compare (typically angle or distance) 
%- neigh: Neighbor proximity (Typically in pixels)
%OUTPUT: 
%-ndist: variable differences within
function [ ndist ] = NearestNbrScat(var, neigh)
ndist = NaN(length(neigh(1,:)),length(neigh(1,:)));

for i = 1: length(neigh(1,:))
    p1 = var(i);    
    for j = 1: length(neigh(1,:))
        p2 = var(j);       
        if j == i
            break
        else             
            ndist(i,j) = abs(p1-p2); 
        end
    end
end

neigh = neigh(:);
ndist = ndist(:);

figure
scatter(neigh,ndist)

end