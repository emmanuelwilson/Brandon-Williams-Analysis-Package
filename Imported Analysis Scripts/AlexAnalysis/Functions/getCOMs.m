function coms = getCOMs(m)
    

    m = m == repmat(nanmax(nanmax(m,[],1),[],2),size(m(:,:,1))); % compute peak location instead

    m = m./repmat(nansum(nansum(m,1),2),[size(m(:,:,1))]);
    
    [mx my] = meshgrid(1:length(m(1,:,1)),1:length(m(:,1,1)));
    x = permute(nansum(nansum(bsxfun(@times,m,mx),1),2),[3 1 2]);
    y = permute(nansum(nansum(bsxfun(@times,m,my),1),2),[3 1 2]);
    coms = [x y];
end