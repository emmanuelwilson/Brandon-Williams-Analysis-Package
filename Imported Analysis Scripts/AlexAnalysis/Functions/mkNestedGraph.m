function mkNestedGraph(d)

    barCount = length(d(:,1)).*nanmax(cellfun(@size,d(:),repmat({2},[numel(d) 1])));

    for di = 1:length(d(:,1))
        for dj = 1:length(d(1,:))
            m = nanmean(d{di,dj},1);
            se = nanstd(d{di,dj},[],1)./sqrt(~isnan(d{di,dj},1));
            patch
        end
    end
end