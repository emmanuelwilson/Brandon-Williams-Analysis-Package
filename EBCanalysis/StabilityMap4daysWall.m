function [StabMap,wallInd] = StabilityMap4daysWall(IndMap,fullpass1W,fullpass2W,fullpass3W,fullpass4W,mouse)

bmap = IndMap;
bmap(bmap>0) = 1;
b = sum(bmap,2);

StabMap = zeros(size(IndMap));
for i = 1 : length(b)
    if ~isempty(find(IndMap(i,1) == fullpass1W,1))
        StabMap(IndMap(i,1),1) = StabMap(IndMap(i,1),1) + 1;
    end
    if ~isempty(find(IndMap(i,2) == fullpass2W,1))
        StabMap(IndMap(i,2),2) = StabMap(IndMap(i,2),2) + 1;
    end
    if ~isempty(find(IndMap(i,3) == fullpass3W,1))
        StabMap(IndMap(i,3),3) = StabMap(IndMap(i,3),3) + 1;
    end
    if ~isempty(find(IndMap(i,4) == fullpass4W,1))
        StabMap(IndMap(i,4),4) = StabMap(IndMap(i,4),4) + 1;
    end    
end
[wallInd(:,1),wallInd(:,2)] = find(StabMap ==1);
figure
imagesc(StabMap)
caxis([0 3])
colorbar
title(mouse)
% legend('No Tuning','Wall Tuning','Object Tuning','Wall and Object Tuning')

end