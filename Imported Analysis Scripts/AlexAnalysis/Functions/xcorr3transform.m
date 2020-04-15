function [ivals tm1 tm2] = xcorr3transform(m1,m2,vals)
    ivals = [];
    m2 = imrotate(m2,vals(3),'nearest','crop');
    
    tm1 = m2(nanmax(vals(1)+1,1):nanmin(length(m1(:,1,1)),length(m1(:,1,1))+vals(1)),...
        nanmax(vals(2)+1,1):nanmin(length(m1(1,:,1)),length(m1(1,:,1))+vals(2)),:);
    tm2 = m1(nanmax(-vals(1)+1,1):nanmin(length(m2(:,1,1)),length(m2(:,1,1))-vals(1)),...
        nanmax(-vals(2)+1,1):nanmin(length(m2(1,:,1)),length(m2(1,:,1))-vals(2)),:);
    
    ivals = nan(length(m1(1,1,:)),1);
    for k = 1:length(m1(1,1,:))
        isGood = ~isnan(tm1(:,:,k))&~isnan(tm2(:,:,k));
        if ~any(isGood(:))
            continue
        end
        km1 = tm1(:,:,k);
        km2 = tm2(:,:,k);
        ivals(k) = corr(km1(isGood),km2(isGood));
    end
%     corr(tm1(repmat(isGood,size(tm1(1,1,:)))),tm2(repmat(isGood,size(tm1(1,1,:)))))
end