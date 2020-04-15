function [pvals actual_mrl null actual_m] = getUnmatchedEBCNMasks(P,HD,T,masks)
    pvals = nan;
    actual = nan;
    null = nan;
    velThresh = 2;
    nsims = 100; %500
    chunkSize = 1;
    minShift = 900;
    
    P = interpNaNs(P')';
    
    [m prepStack] = mkEBCMap2(P,HD,T,masks);
    actual_m = m;
    
    actual_mrl = ebcScore(m);

    null = nan([length(actual_mrl) nsims]);
    parfor iter = 1:nsims./chunkSize
        newGT = repmat(T,[chunkSize 1 1]);
        newGT = circshift(newGT,[0 minShift+randi(length(T(1,:))-minShift.*2)]);
        
        m = mkEBCMap2(P,HD,newGT,masks,[],prepStack);
        
        null(:,iter) = ebcScore(m);
    end
    
    pvals = 1-nanmean(bsxfun(@gt,actual_mrl,null),2);
    pvals(isnan(actual_mrl)) = nan;
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