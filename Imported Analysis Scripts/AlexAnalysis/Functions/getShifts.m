function [peak_fs com_fs xc_fs] = getShifts(m1,m2)
    coms_fs = nan;
    xc_fs = nan;
    peak_fs = nan;
    
    [x y] = meshgrid(1:length(m1(:,1,1)),1:length(m2(1,:,1)));
    nm1 = m1./repmat(nansum(nansum(m1,1),2),size(m1(:,:,1)));
    nm2 = m2./repmat(nansum(nansum(m2,1),2),size(m2(:,:,1)));
    com_m1 = [nansum(nansum(nm1.*repmat(x,size(nm1(1,1,:))),1),2) ...
        nansum(nansum(nm1.*repmat(y,size(nm1(1,1,:))),1),2)];
    com_m2 = [nansum(nansum(nm2.*repmat(x,size(nm1(1,1,:))),1),2) ...
        nansum(nansum(nm2.*repmat(y,size(nm1(1,1,:))),1),2)];
    com_fs = permute(sqrt(nansum([com_m1-com_m2].^2,2)),[3 1 2]);
    
    pm1 = m1==repmat(nanmax(nanmax(m1,[],1),[],2),size(m1(:,:,1)));
    pm2 = m2==repmat(nanmax(nanmax(m2,[],1),[],2),size(m2(:,:,1)));
    peak_fs = nan(length(m1(1,1,:)),1);
    for k = 1:length(m1(1,1,:))
        [x1 y1] = find(pm1(:,:,k));
        [x2 y2] = find(pm2(:,:,k));
        peak_fs(k) = sqrt(nansum([nanmedian(x1)-nanmedian(x2) ...
            nanmedian(y1)-nanmedian(y2)].^2,2));
    end
end