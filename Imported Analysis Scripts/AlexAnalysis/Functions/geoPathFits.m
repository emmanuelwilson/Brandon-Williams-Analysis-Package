function vals = geoPathFits(p1,p2,gT1,gT2)
    vals = [];
    [p1x p2x p1y p2y] = fitDeformationPaths(p1,p2);
    
    m(:,:,:,1) = mkTraceMaps(p1x,gT1,[],[17 17]);
    m(:,:,:,4) = mkTraceMaps(p2x,gT2,[],[17 17]);
    m(:,:,:,2) = mkTraceMaps(p1y,gT1,[],[17 17]);
    m(:,:,:,5) = mkTraceMaps(p2y,gT2,[],[17 17]);
    m(:,:,:,3) = nanmean(m(:,:,:,[1 2]),4);
    m(:,:,:,6) = nanmean(m(:,:,:,[4 5]),4);
    
    vals = nan([3 3 length(m(1,1,:,1))]);
    for i = 1:3
        for j = 4:6
            tm1 = m(:,:,:,i);
            tm2 = m(:,:,:,j);
            rm1 = reshape(tm1,[prod(size(m(:,:,1,1))) length(m(1,1,:,1))]);
            rm2 = reshape(tm2,[prod(size(m(:,:,1,1))) length(m(1,1,:,1))]);
            isBad = any(isnan(rm1),2) | any(isnan(rm2),2);
            tv = corr(rm1(~isBad,:),rm2(~isBad,:));
            vals(i,j-3,:) = tv(logical(eye(size(tv))));
        end
    end
    vals = permute(nanmax(nanmax(vals,[],1),[],2),[3 1 2]);
end