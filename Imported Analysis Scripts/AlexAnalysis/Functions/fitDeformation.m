function [ivals bestFit] = fitDeformation(m1,m2)
    ivals = nan;
    
    sm1 = ~isnan(m1(:,:,1));
    sm2 = ~isnan(m2(:,:,1));
    
    [x1 y1] = meshgrid(1:length(m1(1,:,1)),1:length(m1(:,1,1)));
    [x2 y2] = meshgrid(1:length(m2(1,:,1)),1:length(m2(:,1,1)));
    
    tmp1 = [sm1.*x1];
    tmp1(tmp1==0) = nan;
    sx1 = [nanmin(tmp1,[],2) nanmax(tmp1,[],2)];
    
    tmp1 = [sm1.*y1];
    tmp1(tmp1==0) = nan;
    sy1 = [nanmin(tmp1,[],1)' nanmax(tmp1,[],1)'];
    
    tmp1 = [sm2.*x2];
    tmp1(tmp1==0) = nan;
    sx2 = [nanmin(tmp1,[],2) nanmax(tmp1,[],2)];
    
    tmp1 = [sm2.*y2];
    tmp1(tmp1==0) = nan;
    sy2 = [nanmin(tmp1,[],1)' nanmax(tmp1,[],1)'];
    
    warpX = [nanmax(sx1(:,1),sx2(:,1)) nanmin(sx1(:,2),sx2(:,2))];
    warpY = [nanmax(sy1(:,1),sy2(:,1)) nanmin(sy1(:,2),sy2(:,2))];
    
    warpX = sort(warpX,2);
    warpY = sort(warpY,2);
    
    nmX1 = nan(size(m1));
    nmX2 = nan(size(m2));
    nmY1 = nan(size(m1));
    nmY2 = nan(size(m2));
    power = 2;
    for i = 1:length(warpX(:,1))
        if isnan(warpX(i,1))
            continue
        end
        
        try
            nmX1(i,warpX(i,1):warpX(i,2),:) = interp1(...
                [[([1:nansum(sm1(i,:))]-1)./nansum(sm1(i,:))].*(warpX(i,2)-warpX(i,1))]',...
                permute(m1(i,sm1(i,:),:),[3 2 1])',[0:warpX(i,2)-warpX(i,1)]','pchip','extrap');
        end
        try
            nmX2(i,warpX(i,1):warpX(i,2),:) = interp1(...
                [[([1:nansum(sm2(i,:))]-1)./nansum(sm2(i,:))].*(warpX(i,2)-warpX(i,1))]',...
                permute(m2(i,sm2(i,:),:),[3 2 1])',[0:warpX(i,2)-warpX(i,1)]','pchip','extrap');   
        end
    end
    
    for i = 1:length(warpY(:,1))
        if isnan(warpY(i,1))
            continue
        end
        try
            nmY1(warpY(i,1):warpY(i,2),i,:) = interp1(...
                [[([1:nansum(sm1(:,i))]-1)./nansum(sm1(:,i))].*(warpY(i,2)-warpY(i,1))]',...
                permute(m1(sm1(:,i),i,:),[3 1 2])',[0:warpY(i,2)-warpY(i,1)]','pchip','extrap');
        end
        try
            nmY2(warpY(i,1):warpY(i,2),i,:) = interp1(...
                [[([1:nansum(sm2(:,i))]-1)./nansum(sm2(:,i))].*(warpY(i,2)-warpY(i,1))]',...
                permute(m2(sm2(:,i),i,:),[3 1 2])',[0:warpY(i,2)-warpY(i,1)]','pchip','extrap');
        end
             
    end
    
    nmC1 = nanmean(cat(4,nmX1,nmY1),4);
    nmC2 = nanmean(cat(4,nmX2,nmY2),4);
    
    %% Geo examples
    
%     for k = 1:length(m1(1,1,:))
%         figure
%         set(gcf,'position',[50 50 800 400])
%         subplot(2,4,1)
%         imagesc(m1(:,:,k))
%         alpha(double(~isnan(m1(:,:,k))))
%         axis off
%         axis equal
%         subplot(2,4,5)
%         imagesc(m2(:,:,k))
%         alpha(double(~isnan(m2(:,:,k))))
%         axis off
%         axis equal
%         subplot(2,4,2)
%         imagesc(nmX1(:,:,k))
%         alpha(double(~isnan(nmX1(:,:,k))))
%         axis off
%         axis equal
%         subplot(2,4,6)
%         imagesc(nmX2(:,:,k))
%         alpha(double(~isnan(nmX2(:,:,k))))
%         axis off
%         axis equal
%         subplot(2,4,3)
%         imagesc(nmY1(:,:,k))
%         alpha(double(~isnan(nmY1(:,:,k))))
%         axis off
%         axis equal
%         subplot(2,4,7)
%         imagesc(nmY2(:,:,k))
%         alpha(double(~isnan(nmY2(:,:,k))))
%         axis off
%         axis equal
%         subplot(2,4,4)
%         imagesc(nmC1(:,:,k))
%         alpha(double(~isnan(nmC1(:,:,k))))
%         axis off
%         axis equal
%         subplot(2,4,8)
%         imagesc(nmC2(:,:,k))
%         alpha(double(~isnan(nmC2(:,:,k))))
%         axis off
%         axis equal
%         colormap jet
%         saveFig(gcf,['Plots/GeoFitExamples/Cell_' num2str(k)],[{'tiff'} {'pdf'}])
%         close all
%     end
    %%
    
    isGood = repmat(~isnan(nmX1(:,:,1)) & ~isnan(nmY1(:,:,1)) & ...
        ~isnan(nmX2(:,:,1)) & ~isnan(nmY2(:,:,1)),[1 1 length(nmX1(1,1,:))]);
    
    allCompMaps = cat(4,nmX1,nmY1,nmC1,nmX2,nmY2,nmC2);
    vals = nan([length(allCompMaps(1,1,1,:)) length(allCompMaps(1,1,1,:)) length(m1(1,1,:))]);
    bestFits = nan([length(allCompMaps(1,1,1,:)) length(allCompMaps(1,1,1,:)) length(m1(1,1,:))]);
    for i = 1:length(allCompMaps(1,1,1,:))
        a = allCompMaps(:,:,:,i);
        for j = i+1:length(allCompMaps(1,1,1,:))
            if i > 3 || j < 4
                continue
            end
            b = allCompMaps(:,:,:,j);
            tmp = corr(reshape(a(isGood),[nansum(nansum(isGood(:,:,1))) length(a(1,1,:))]),...
                reshape(b(isGood),[nansum(nansum(isGood(:,:,1))) length(b(1,1,:))]));
            vals(i,j,:) = tmp(logical(eye(size(tmp))));
        end
    end
    ivals = permute(nanmax(nanmax(vals,[],1),[],2),[3 1 2]);
    
    tmp = vals(1:3,4:6,:);
    inds = find(tmp==repmat(nanmax(nanmax(tmp,[],1),[],2),[size(tmp(:,:,1))]));
    [a b c] = ind2sub(size(tmp),inds);
    
    area1 = nansum(nansum(~isnan(nanmax(m1,[],3))));
    area2 = nansum(nansum(~isnan(nanmax(m2,[],3))));
    if area1 > area2
        bestFit = [a];
    else
        bestFit = [b];
    end
end




















