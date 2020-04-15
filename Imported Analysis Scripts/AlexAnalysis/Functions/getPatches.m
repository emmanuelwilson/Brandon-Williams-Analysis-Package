function [labelMap patchd2c patchCenter patchMax] = getPatches(map,minSize,a,b)
    labelMap = zeros(size(map));
    fieldNum = 0;
    minAct = a;
    if isempty(minAct)
        minAct = (b.*nanmax(map(:)));
    end
    while sum(labelMap(:)==0 & map(:)>=minAct)>0
        fieldNum = fieldNum+1;
        [cx cy] = find(labelMap==0 & map>=minAct,1,'first');
        while ~isempty([cx cy])
            labelMap(sub2ind(size(map),cx,cy)) = fieldNum;
            [cx cy] = find(labelMap==fieldNum);
            poss = [[cx+1 cy]; [cx cy+1]; [cx-1 cy]; [cx cy-1]];
            poss = poss((poss(:,1)<=length(map(:,1))&poss(:,1)>0 & poss(:,2)<=length(map(1,:))&poss(:,2)>0),:);
            poss = poss(map(sub2ind(size(map),poss(:,1),poss(:,2)))>minAct & ...
                labelMap(sub2ind(size(map),poss(:,1),poss(:,2)))==0,:);
            cx = poss(:,1);
            cy = poss(:,2);
        end
    end
    isField = false(size(map));
    for i = 1:fieldNum
        if sum(labelMap(:)==i)>=minSize
            isField(labelMap==i) = true;
        end
    end
    
    labelMap(~isField)=0;
    iter = 0;
    for i = unique(labelMap)'
        labelMap(labelMap==i)=iter;
        iter = iter+1;
    end
    
    % Order by distance to center
    [x y] = meshgrid(1:length(labelMap(1,:)),1:length(labelMap(:,1)));
    x = x-ceil(length(labelMap(1,:))./2);
    y = y-ceil(length(labelMap(:,1))./2);
    d2c = sqrt(x.^2 + y.^2);
    patchDists = [];
    patchCenter = [];
    patchMax = [];
    for i = unique(labelMap(labelMap>0))'
        
        tmp = (labelMap==i);
        tmp = tmp./sum(tmp(:));
        [x2 y2] = meshgrid(1:length(tmp(1,:)),1:length(tmp(:,1)));
        com = [sum(tmp(:).*x2(:)) sum(tmp(:).*y2(:))];
        patchCenter(i,:) = com;
        
        tmp = (labelMap==i);
        [x y] = find(map == nanmax(map(tmp)) & tmp);
        patchMax(i,:) = [ceil(median(x)) ceil(median(y))];
        
        tmp = map;
        tmp(labelMap~=i) = 0;
        tmp = tmp./sum(tmp(:));
        tmp = tmp.*d2c;
        patchDists(i) = sum(tmp(:));
        
        patchDists(i) = d2c(patchMax(i,1),patchMax(i,2));
    end
    [a b] = sort(patchDists);
    newLabelMap = labelMap;
    for i = unique(labelMap(labelMap>0))'
        newLabelMap(labelMap==i) = find(b==i);
    end
    labelMap = newLabelMap;
    patchd2c = a;
    patchCenter = patchCenter(b,:);
    patchMax = patchMax(b,:);
end