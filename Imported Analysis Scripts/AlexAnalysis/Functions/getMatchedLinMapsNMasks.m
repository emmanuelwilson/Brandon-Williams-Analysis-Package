function [m s val ival] = getMatchedLinMapsNMasks(p,ts,masks)
    ival = nan(length(masks));
    val = nan(length(masks));
    for k = 1:length(masks)
        [m(:,:,:,k) s(:,:,k)] = mkLinMaps(p,ts,masks{k});
    end
    
    os = min(s,[],3);
    
%     ival = nan([length(masks) length(masks) length(m(1,1,:,1))]);
%     val = nan(length(masks));
%     for ki = 1:length(masks)
%         for kj = ki+1:length(masks)
%             ma = m(:,:,:,ki);
%             mb = m(:,:,:,kj);
%             val(ki,kj) = corr(ma(~isnan(ma)&~isnan(mb)),...
%                 mb(~isnan(ma)&~isnan(mb)));
%             if nargout > 3
%                 
%                 mak = reshape(ma,numel(ma(:,:,1)),[]);
%                 mbk = reshape(mb,numel(mb(:,:,1)),[]);
%                 tmp = corr(mak(any(~isnan(mak),2)&any(~isnan(mbk),2),:),...
%                     mbk(any(~isnan(mak),2)&any(~isnan(mbk),2),:));
%                 ival(ki,kj,:) = (tmp(logical(eye(length(tmp)))));
%             end
%         end
%     end 

    
    
    nsim = 100;
    if all(os(:)==0)
        val = nan;
    	return
    end
    
    for k = 1:length(masks)
        om(:,:,:,:,k) = mkLinMaps(p,ts,masks{k},[],os,nsim);
    end
    
    ival = nan([length(masks) length(masks) length(ts(:,1))]);
    
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