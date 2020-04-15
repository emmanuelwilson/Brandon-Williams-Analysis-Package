function [m os val ival] = getMatchedMapsNMasks(p,ts,masks)
    ival = nan;
    
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
        om(:,:,:,:,k) = mkTraceMaps(p,ts,masks{k},[],os,nsim);
    end
        
    ival = nan([length(masks) length(masks) length(ts(:,1))]);
    
    for ki = 1:length(masks)
        for kj = ki+1:length(masks)
            ma = om(:,:,:,:,ki);
            mb = om(:,:,:,:,kj);

            a = reshape(ma,[numel(ma(:,:,:,1)) nsim]);
            b = reshape(mb,[numel(ma(:,:,:,1)) nsim]);
            a(isnan(a(:,1)),:) = [];
            b(isnan(b(:,1)),:) = [];
            
            tmp = corr(a,b);
            val(ki,kj) = nanmedian(tmp(logical(eye(length(tmp)))));
                        
            if nargout > 3
                for k = 1:length(ma(1,1,:,1))
                    a = reshape(ma(:,:,k,:),[numel(ma(:,:,1,1)) nsim]);
                    b = reshape(mb(:,:,k,:),[numel(ma(:,:,1,1)) nsim]);
                    a(isnan(a(:,1)),:) = [];
                    b(isnan(b(:,1)),:) = [];
                    tmp = corr(a,b);
                    ival(ki,kj,k) = nanmean(tmp(logical(eye(length(tmp)))));
                end
            end
        end
    end
end