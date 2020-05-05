function [mrl best_d pfd] = ebcScore(m,wm)
    angF = repmat([2.*pi./(length(m(:,1)).*2):2.*pi./length(m(:,1)):2.*pi]',[1 length(m(1,:,1))]);
    aR = nan(length(m(1,1,:)),length(m(1,:,1)));
%     aWR = nan(length(m(1,1,:)),length(m(1,:,1)));
    
    medianCorr = nan(length(m(1,1,:)),1);
    pos = nan(length(m(1,1,:)),2);
    for k = 1:length(m(1,1,:))
%         tmp = nan(length(wm(1,1,1,:)),length(wm(1,:,1,1)));
%         for j = 1:length(wm(1,1,1,:))
%             tmp(j,:) = circ_r(angF,wm(:,:,k,j),[],1);
%         end
%         aWR(k,:) = nansum(log(tmp),1);
        
        aR(k,:) = circ_r(angF,m(:,:,k),[],1);
        [x y] = find(m(:,:,k)==nanmax(nanmax(m(:,:,k))),1,'first');
        pos(k,:) = [x y];
        if nargin > 1 && ~isempty(wm)
            tmp = reshape(wm(:,:,k,:),[numel(wm(:,:,1,1)) length(wm(1,1,1,:))]);
            tmp = corr(tmp(~any(isnan(tmp),2),:),tmp(~any(isnan(tmp),2),:));
            medianCorr(k) = nanmedian(tmp(logical(triu(true(size(tmp)),1))));
        end
    end
    
    [mrl b] = nanmax(aR,[],2);
    
    [a b] = nanmax(permute(nanmean(m,1),[3 2 1]),[],2);
    
%     mrl = nan(length(m(1,1,:)),1);
%     best_d = b;
%     for i = 1:length(b)
%         mrl(i) = aR(i,b(i));
%     end
    
    best_d = pos(:,2);
    pfd = (pos(:,1));
    
%     tmp = permute(nanmean(m,1),[3 2 1]);
%     tmp = bsxfun(@minus,tmp,nanmin(tmp,[],2));
%     tmp = bsxfun(@rdivide,tmp,nansum(tmp,2));
%     nansum(bsxfun(@times,tmp,1:length(tmp(1,:))),2)
    
%     tm = m-repmat(nanmin(nanmin(m,[],1),[],2),size(m(:,:,1)));
%     nm = tm./repmat(nansum(nansum(tm,1),2),size(tm(:,:,1)));
%     [x y] = meshgrid(1:length(m(1,:,1)),1:length(m(:,1,1)));
%     
%     nansum(nansum(bsxfun(@times,nm,x),1),2)
end