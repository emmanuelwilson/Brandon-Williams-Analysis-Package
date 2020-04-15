function [pvals actual_shc null actual_m] = getshEBC(P,HD,T)
    pvals = nan;
    actual = nan;
    null = nan;
    velThresh = 2;
    nsims = 100; %500
    chunkSize = 1;
    minShift = 900;
    
    P = interpNaNs(P')';

    [actual_shc actual_m] = help_sh(P,HD,T);
    
    
    null = nan([length(actual_shc) nsims]);
    parfor iter = 1:nsims./chunkSize
        newGT = repmat(T,[chunkSize 1 1]);
        newGT = circshift(newGT,[0 minShift+randi(length(T(1,:))-minShift.*2)]);
        
        [null(:,iter) m] = help_sh(P,HD,newGT);
    end
    
    pvals = 1-nanmean(bsxfun(@gt,actual_shc,null),2);
    pvals(isnan(actual_shc)) = nan;
end

function [xc m] = help_sh(P,HD,T)
    mask = ([1:length(HD)]<length(HD)./2);
    m1 = mkEBCMap2(P,HD,T,{mask});
    m2 = mkEBCMap2(P,HD,T,{~mask});

    rm1 = reshape(m1,[numel(m1(:,:,1)) length(m1(1,1,:))]);
    rm2 = reshape(m2,[numel(m2(:,:,1)) length(m2(1,1,:))]);
    isGood = ~any(isnan(rm1),2) & ~any(isnan(rm2),2);
    xc = corr(rm1(isGood,:),rm2(isGood,:));
    xc = xc(logical(eye(size(xc))));
    
    m = cat(4,m1,m2);
end