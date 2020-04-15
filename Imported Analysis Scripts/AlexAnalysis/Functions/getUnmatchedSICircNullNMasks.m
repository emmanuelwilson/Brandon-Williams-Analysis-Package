function [pvals actual null ] = getUnmatchedSICircNullNMasks(P,T,masks)
    pvals = nan;
    actual = nan;
    null = nan;
    velThresh = 2;
    nsims = 1000; %500
    minShift = 900;
    
    actual = help_getSI(P,T,masks);
    
    null = nan([length(actual) nsims]);
    parfor si = 1:nsims
        gT = circshift(T,[0 minShift+randi(length(T(1,:))-minShift.*2)]);
        null(:,si) = help_getSI(P,gT,masks);
    end
    
    pvals = 1-nanmean(bsxfun(@gt,actual,null),2);
    pvals(isnan(actual)) = nan;
end

function si = help_getSI(P,T,masks)
    for k = 1:length(masks)
        [blah(:,:,:,k) s(:,:,k) m(:,:,:,k)] = mkTraceMaps(P,T,masks{k});
    end
    m = reshape(m,[numel(m(:,:,1)) length(m(1,1,:))]);
    s = s(:);
    m(s<15,:) = [];
    s(s<15) = [];
    normS = s./nansum(s(:));
    normM = bsxfun(@rdivide,m,nanmean(m,1));
    si = nansum(bsxfun(@times,normS,normM.*log(normM)),1)';
end