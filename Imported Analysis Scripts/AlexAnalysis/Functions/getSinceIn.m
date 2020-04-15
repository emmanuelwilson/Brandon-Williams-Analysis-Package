function [isMostRecent indexSinceIn distanceSinceIn] = getSinceIn(isIn,p)
    isMostRecent = isIn;
    for i = 2:length(isIn(1,:))
        if all(isIn(:,i)==0)
            isMostRecent(:,i) = isMostRecent(:,i-1);
        end
    end
    
    inBetween = any(isMostRecent&~isIn,1);
    distanceSinceIn = nan(1,length(inBetween));
    while any(inBetween)
        start = find(inBetween,1,'first');
        stop = find(~inBetween(start+1:end),1,'first')-1;
        if isempty(stop)
            stop = length(inBetween)-start;
        end
        tmp = nanmax(sqrt(nansum(bsxfun(@minus,p(:,start:start+stop),p(:,start)).^2,1)));
%         distanceSinceIn(start:start+stop) = tmp;
        distanceSinceIn(start+stop) = tmp;
        inBetween(start:start+stop) = false;
    end
    
    indexSinceIn = zeros(size(isIn));
    for j = 1:length(isIn(:,1))
        for i = 2:length(isIn(1,:))
            if all(isIn(j,i)==0)
                indexSinceIn(j,i) = indexSinceIn(j,i-1)+1;
            end
        end
    end
end