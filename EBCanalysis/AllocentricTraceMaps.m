function [m s unm] = AllocentricTraceMaps(p,t,mask,dims,downsample,nsim)
    if nargin<3 || isempty(mask)
        mask = true(length(p(1,:)));
    end
    
    if nargin<5 || isempty(downsample)
        downsample = [];
    end
    if nargin<5 || isempty(downsample) || nargin<6 || isempty(nsim)
        nsim = 1;
    end
    
    
%     t = [nan(length(t(:,1)),1) diff(t,[],2)>0];
%     rising = trimNoise(rising);
    
    binsize = 7; % 2.5
    kern = 10; %    4
    
    p = bsxfun(@minus,p,nanmin(p')');
    p = floor(p./binsize)+1;
    if nargin<4 || isempty(dims)
        m = nan([nanmax(p') length(t(:,1)) nsim]);
    else
        m = nan([dims length(t(:,1)) nsim]);
    end
    s = zeros(size(m(:,:,1)));
    
    p(:,~mask) = [];
    t(:,~mask) = [];
    
    up = unique(p','rows');
    up(any(isnan(up),2),:) = [];
    for loc = up'
        isGood = find([loc(1)==p(1,:) & loc(2)==p(2,:)]);        
        m(loc(1),loc(2),:) = nansum(t(:,isGood),2);
        s(loc(1),loc(2)) = length(isGood);
    end
    unm = m;
    bad = isnan(m);
    m(bad) = 0;    
    m = imfilter(m,fspecial('gauss',round([kern.*5 kern.*5]./binsize),kern./binsize),'same');
    m(bad) = nan;

    
%     m = bsxfun(@times,m,s);
%     m = imfilter(m,fspecial('gauss',round([kern.*5 kern.*5]./binsize),kern./binsize),'same');
%     s = imfilter(s,fspecial('gauss',round([kern.*5 kern.*5]./binsize),kern./binsize),'same');
%     m = bsxfun(@rdivide,m,s);
end