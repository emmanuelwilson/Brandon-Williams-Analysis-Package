function [m sa val ival] = getMatchedSamplingMaps(p,ts,maskA,maskB)
    ival = nan;
    
    [ma sa] = mkTraceMaps(p,ts,maskA);
    [mb sb] = mkTraceMaps(p,ts,maskB);
    
    m = cat(4,ma,mb);
    
    os = min(sa,sb);
    
%     figure(1)
%     set(gcf,'position',[50 50 500 500])
%     subplot(2,2,1)
%     imagesc(sa)
%     colormap hot
%     alpha(double(~isnan(ma(:,:,1))))
%     caxis([0 nanmax([sa(:); sb(:)])])
%     axis off
%     subplot(2,2,2)
%     imagesc(sb)
%     colormap hot
%     alpha(double(~isnan(mb(:,:,1))))
%     caxis([0 nanmax([sa(:); sb(:)])])
%     axis off
%     for i = 3:4
%         subplot(2,2,i)
%         imagesc(os)
%         colormap hot
%         alpha(double(os~=0))
%         caxis([0 nanmax([sa(:); sb(:)])])
%         axis off
%     end
    
    nsim = 100;
    if all(os(:)==0)
        val = nan;
    	return
    end
    [ma] = mkTraceMaps(p,ts,maskA,[],os,nsim);
    [mb] = mkTraceMaps(p,ts,maskB,[],os,nsim);
    a = reshape(ma,[numel(ma(:,:,:,1)) nsim]);
    b = reshape(mb,[numel(ma(:,:,:,1)) nsim]);
    a(isnan(a(:,1)),:) = [];
    b(isnan(b(:,1)),:) = [];
    tmp = corr(a,b);
    val = nanmean(tmp(logical(eye(length(tmp)))));
    
    if nargout > 3
        ival = nan(length(ma(1,1,:,1)),1);
        for k = 1:length(ma(1,1,:,1))
            a = permute(ma(:,:,k,:),[1 2 4 3]);
            b = permute(mb(:,:,k,:),[1 2 4 3]);
            a = reshape(a,[numel(a(:,:,1)) nsim]);
            b = reshape(b,[numel(b(:,:,1)) nsim]);
            a(isnan(a(:,1)),:) = [];
            b(isnan(b(:,1)),:) = [];
            tmp = corr(a,b);
            ival(k) = nanmean(tmp(logical(eye(length(tmp)))));
        end
    end
end