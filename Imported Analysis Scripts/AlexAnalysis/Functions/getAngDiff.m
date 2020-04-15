function [ad a] = getAngDiff(m1,m2)
    axc = nan(length(m1(:,1,1)),length(m1(1,1,:)));
    for i = 1:length(m1(:,1,1))
        tmp1 = m1;
        tmp2 = m2;

        tmp2 = circshift(tmp2,[i-1 0 0]);

        rt1 = reshape(tmp1,[numel(tmp1(:,:,1)) length(tmp1(1,1,:))]);
        rt2 = reshape(tmp2,[numel(tmp2(:,:,1)) length(tmp2(1,1,:))]);

        isGP = ~[any(isnan(rt1),2)|any(isnan(rt2),2)];

        xc = corr(rt1(isGP,:),rt2(isGP,:));
        axc(i,:) = xc(logical(eye(size(xc))));
    end
    [a b] = nanmax(axc,[],1);    
    bestRot = mode(b);
    
%     axc = bsxfun(@rdivide,bsxfun(@minus,axc,nanmin(axc,[],1)),range(axc,1));
%     [blah bestRot] = nanmax(nanmean(axc,2));
    a = axc(bestRot,:)';
    
%     b = b-bestRot;
    ad = mod((b-1).*360./length(axc(:,1)),360);
    ad = nanmin(ad,360-ad)';
end