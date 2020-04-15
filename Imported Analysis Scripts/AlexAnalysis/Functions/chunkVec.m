function cv = chunkVec(v,t)
    cv = nan(length(v(:,1)),ceil(length(v(1,:))./t));
    for i = 1:t:length(v(1,:))
        cv(:,((i-1)./t)+1) = nanmean(v(:,i:nanmin(i+t,length(v(1,:)))),2);
    end
end