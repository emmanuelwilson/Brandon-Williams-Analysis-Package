function [m s] = mkHDMaps(p,t,mask,dims,downsample,nsim)
    if nargin<3 || isempty(mask)
        mask = true(1,length(p(1,:)));
    end
    
    if nargin<5 || isempty(downsample)
        downsample = [];
    end
    if nargin<5 || isempty(downsample) || nargin<6 || isempty(nsim)
        nsim = 1;
    end
    
    
%     t = [nan(length(t(:,1)),1) diff(t,[],2)>0];
%     rising = trimNoise(rising);
    
    binsize = 8;
    kern = 12; %3
    
    np = mod(180.*p./(pi),360);
    inp = interpNaNs(np');
    p = floor(inp./binsize)+1;

    p(~mask) = [];
    t(:,~mask) = [];

    up = unique(p);
    up(any(isnan(up),2),:) = [];
    
    m = nan(nanmax(up),length(t(:,1)));
    s = nan(nanmax(up),1);
    for loc = up'
        isGood = find([loc==p]);
        if ~isempty(downsample)
            repGood = nan(nsim,downsample(loc(1),loc(2)));
            for si = 1:nsim
                repGood(si,:) = isGood(randperm(length(isGood),downsample(loc(1),loc(2))));
            end
            repGood = repGood';
            tmp = reshape(t(:,repGood),[length(t(:,1)) downsample(loc(1),loc(2)) nsim]);
            m(loc(1),loc(2),:,:) = permute(nanmean(tmp,2),[4 2 1 3]);
        else
            m(loc,:) = nanmean(t(:,isGood),2);
        end
        s(loc) = length(isGood);
    end
    
    m = repmat(m,[3 1]);
    bad = isnan(m);
    m(bad) = 0;
    m = imfilter(m,fspecial('gauss',round([kern.*5 binsize]./binsize)),'same');
    m(bad) = nan;
    m = m(length(s)+1:2.*length(s),:);
end