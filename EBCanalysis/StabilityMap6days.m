function [StabMap,wallInd,obInd,wallobInd] = StabilityMap6days(IndMap,fullpass1O,fullpass2O,fullpass3O,fullpass4O,fullpass5O,fullpass6O,fullpass1W,fullpass2W,fullpass3W,fullpass4W,fullpass5W,fullpass6W,mouse)

bmap = IndMap;
bmap(bmap>0) = 1;
b = sum(bmap,2);

StabMap = zeros(size(IndMap));
for i = 1 : length(b)
    %     if b(i) == 4
    if ~isempty(find(IndMap(i,1) == fullpass1W,1))
        StabMap(IndMap(i,1),1) = StabMap(IndMap(i,1),1) + 1;
    end
    if ~isempty(find(IndMap(i,1) == fullpass1O,1))
        StabMap(IndMap(i,1),1) = StabMap(IndMap(i,1),1) + 2;
    end
    if ~isempty(find(IndMap(i,2) == fullpass2W,1))
        StabMap(IndMap(i,2),2) = StabMap(IndMap(i,2),2) + 1;
    end
    if ~isempty(find(IndMap(i,2) == fullpass2O,1))
        StabMap(IndMap(i,2),2) = StabMap(IndMap(i,2),2) + 2;
    end
    if ~isempty(find(IndMap(i,3) == fullpass3W,1))
        StabMap(IndMap(i,3),3) = StabMap(IndMap(i,3),3) + 1;
    end
    if ~isempty(find(IndMap(i,3) == fullpass3O,1))
        StabMap(IndMap(i,3),3) = StabMap(IndMap(i,3),3) + 2;
    end
    if ~isempty(find(IndMap(i,4) == fullpass4W,1))
        StabMap(IndMap(i,4),4) = StabMap(IndMap(i,4),4) + 1;
    end
    if ~isempty(find(IndMap(i,4) == fullpass4O,1))
        StabMap(IndMap(i,4),4) = StabMap(IndMap(i,4),4) + 2;
    end
    if ~isempty(find(IndMap(i,5) == fullpass5W,1))
        StabMap(IndMap(i,5),5) = StabMap(IndMap(i,5),5) + 1;
    end
    if ~isempty(find(IndMap(i,5) == fullpass5O,1))
        StabMap(IndMap(i,5),5) = StabMap(IndMap(i,5),5) + 2;
    end
    if ~isempty(find(IndMap(i,6) == fullpass6W,1))
        StabMap(IndMap(i,6),6) = StabMap(IndMap(i,6),6) + 1;
    end
    if ~isempty(find(IndMap(i,6) == fullpass6O,1))
        StabMap(IndMap(i,6),6) = StabMap(IndMap(i,6),6) + 2;
    end
%     end
end
[wallInd(:,1),wallInd(:,2)] = find(StabMap ==1);
[obInd(:,1),obInd(:,2)] = find(StabMap ==2);
[wallobInd(:,1),wallobInd(:,2)] = find(StabMap ==3);
figure
imagesc(StabMap)
caxis([0 3])
colorbar
title(mouse)
% legend('No Tuning','Wall Tuning','Object Tuning','Wall and Object Tuning')

end