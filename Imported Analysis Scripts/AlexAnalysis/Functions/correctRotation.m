function m = correctRotation(m)
    om = m;
    for j = 2:length(m(1,1,1,:))
        tmp1 = om(:,:,:,1);
        tmp2 = om(:,:,:,j);
        vals = nan(1,4);
        for rot = 0:3
            rtmp2 = imrotate(tmp2,rot.*90);
            goodPixels = ~isnan(rtmp2(:,:,1))&~isnan(tmp1(:,:,1));
            vals(rot+1) = corr(tmp1(repmat(goodPixels,[1 1 length(tmp1(1,1,:))])), ...
                rtmp2(repmat(goodPixels,[1 1 length(tmp1(1,1,:))])));
        end
        [a bestRot] = nanmax(vals);
        m(:,:,:,j) = imrotate(tmp2,(bestRot-1).*90);
    end
end