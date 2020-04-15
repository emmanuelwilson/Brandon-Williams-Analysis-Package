function [pvals actual null ] = getUnmatchedCircNullNMasks(P,T,masks)
    pvals = nan;
    actual = nan;
    null = nan;
    velThresh = 2;
    nsims = 1000; %500
    minShift = 900;
    
    actual = help_getMaskCorrs(P,T,masks);
    
    null = nan([size(actual) nsims]);
    parfor si = 1:nsims
        gT = circshift(T,[0 minShift+randi(length(T(1,:))-minShift.*2)]);
        null(:,:,:,si) = help_getMaskCorrs(P,gT,masks);
    end
    
    pvals = 1-nanmean(bsxfun(@gt,actual,null),4);
    pvals(isnan(actual)) = nan;
end

function ival = help_getMaskCorrs(P,T,masks)
    for k = 1:length(masks)
        [m(:,:,:,k) s(:,:,k)] = mkTraceMaps(P,T,masks{k});
    end
    
    ival = nan([length(masks) length(masks) length(T(:,1))]);
    for ki = 1:length(masks)
        for kj = ki+1:length(masks)
            ma = m(:,:,:,ki);
            mb = m(:,:,:,kj);
            
            a = reshape(ma(:,:,:),[numel(ma(:,:,1,1)) length(ma(1,1,:,1))]);
            b = reshape(mb(:,:,:),[numel(ma(:,:,1,1)) length(ma(1,1,:,1))]);
            isGood = ~isnan(a(:,1))&~isnan(b(:,1));
            a(~isGood,:) = [];
            b(~isGood,:) = [];
            
            tmp = corr(a,b);
            ival(ki,kj,:) = tmp(logical(eye(length(tmp))));
        end
    end
end