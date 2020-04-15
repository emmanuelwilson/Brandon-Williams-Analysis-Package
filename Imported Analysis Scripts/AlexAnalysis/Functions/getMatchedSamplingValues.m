function [val ival mfr pfr] = getMatchedSamplingValues(p,ts,masks)
    ival = nan;
    mfr = nan;
    pfr = [];
    pfshift = [];
    for k = 1:length(masks)
        [m(:,:,:,k) s(:,:,k)] = mkTraceMaps(p,ts,masks{k});
    end
    
    os = min(s,[],3);
    
    
    nsim = 100;
    if all(os(:)==0)
        val = nan;
    	return
    end
    
    for k = 1:length(masks)
        [om(:,:,:,:,k) uns(:,:,:,:,k) unm(:,:,:,:,k)] = mkTraceMaps(p,ts,masks{k},[],os,nsim);
    end
    
%     unm = om;
    unm(:,:,:,:,1:2) = bsxfun(@rdivide,unm(:,:,:,:,1:2),nanmax(nansum(nansum(unm(:,:,:,:,1:2),1),2),[],5)); % normalize
    unm(:,:,:,:,3:4) = bsxfun(@rdivide,unm(:,:,:,:,3:4),nanmax(nansum(nansum(unm(:,:,:,:,3:4),1),2),[],5)); % normalize

%     unm = bsxfun(@rdivide,unm,nanmax(nanmax(nanmax(unm,[],1),[],2),[],5)); % normalize
    
%     ros = repmat(os,[1 1 length(om(1,1,:,1,1)) nsim]);

    ival = nan([length(masks) length(masks) length(ts(:,1))]);
    mfr = nan([length(masks) length(masks) length(ts(:,1))]);
    pfr = nan([length(masks) length(masks) length(ts(:,1))]);
    pfr_shift = nan([length(masks) length(masks) length(ts(:,1))]);
    for ki = 1:length(masks)
        for kj = ki+1:length(masks)
            ma = om(:,:,:,:,ki);
            mb = om(:,:,:,:,kj);
            
            %%% norm to max to eliminate relative weight
%             normA = nanmax(nanmax(ma,[],1),[],2);
%             normA(normA==0) = 1;
%             normB = nanmax(nanmax(mb,[],1),[],2);
%             normB(normB==0) = 1;
% 
%             ma = bsxfun(@rdivide,ma,normA);
%             mb = bsxfun(@rdivide,mb,normB);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            a = reshape(ma,[numel(ma(:,:,:,1)) nsim]);
            b = reshape(mb,[numel(ma(:,:,:,1)) nsim]);
            a(isnan(a(:,1)),:) = [];
            b(isnan(b(:,1)),:) = [];
            tmp = corr(a,b);
            val(ki,kj) = nanmean(tmp(logical(eye(length(tmp)))));
            
            if nargout > 1
                for k = 1:length(ma(1,1,:,1))
                    a = reshape(ma(:,:,k,:),[numel(ma(:,:,1,1)) nsim]);
                    b = reshape(mb(:,:,k,:),[numel(ma(:,:,1,1)) nsim]);
                    a(isnan(a(:,1)),:) = [];
                    b(isnan(b(:,1)),:) = [];
                    tmp = corr(a,b);
                    ival(ki,kj,k) = nanmean(tmp(logical(eye(length(tmp)))));
                end
            end
            
            if nargout > 2
                sa = nansum(nansum(unm(:,:,:,:,ki),1),2);
                sb = nansum(nansum(unm(:,:,:,:,kj),1),2);
                mfr(ki,kj,:) = nanmean(abs(sa-sb),4);
            end
            
            if nargout > 3
                sa = nanmax(nanmax(unm(:,:,:,:,ki),[],1),[],2);
                sb = nanmax(nanmax(unm(:,:,:,:,kj),[],1),[],2);
                pfr(ki,kj,:) = nanmean(abs(sa-sb),4);
            end
        end
    end
end