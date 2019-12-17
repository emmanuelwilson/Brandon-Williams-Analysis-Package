%Putting together truth map
addon = zeros(length(truthmap(:,1)),1);
truthmap = cat(2,truthmap,addon);
temp = cell_registered_struct6.cell_to_index_map;
for i = 1 :485
    loc = find(temp(:,1) == i);
    if temp(loc,2) > 0
        loct = find(truthmap(:,1) == i);
        truthmap(loct,end) = temp(loc,2);
    end
end