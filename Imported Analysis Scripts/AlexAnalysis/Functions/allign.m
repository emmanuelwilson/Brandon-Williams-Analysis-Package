function ts = allign(ts1,ts2,lag)
    if ~exist('lag','var')
        lag = 0;
    end
    ts2 = ts2+lag;
    ts = nan(1,length(ts2));
    dists = nan(length(ts2),1);
    for t = ts2
        [dists(ts2==t) ts(ts2==t)] = min(abs(ts1-t));
    end
    ts(dists>median(diff(ts1))) = [];
end