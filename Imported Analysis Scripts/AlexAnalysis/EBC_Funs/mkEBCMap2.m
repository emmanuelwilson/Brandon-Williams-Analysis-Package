function [m prepStack] = mkEBCMap2(P,HD,T,mask,walls,prepStack)
    
    params = ebcMapParams();

    kern = params.kern;
    binDist = params.binDist;
    binAng = params.binAng;
    stepSize = params.stepSize;
    
    if nargin >=6 && ~isempty(prepStack)
        
        d = prepStack.d;
        t = prepStack.t;
        uLoc = prepStack.uLoc;
        walls = prepStack.walls;
    else
        wm = [];
        if nargin < 5 || isempty(walls)
            walls = [0 0 0 nanmax(P(:)); 0 0 nanmax(P(:)) 0; ...
                    0 nanmax(P(:)) nanmax(P(:)) nanmax(P(:)); ...
                    nanmax(P(:)) 0 nanmax(P(:)) nanmax(P(:))];

            walls = cat(3,[zeros(1,[nanmax(P(:))./stepSize]+1); 0:stepSize:nanmax(P(:))],...
                [0:stepSize:nanmax(P(:)); zeros(1,[nanmax(P(:))./stepSize]+1)], ...
                [ones(1,[nanmax(P(:))./stepSize]+1).*nanmax(P(:)); 0:stepSize:nanmax(P(:))] ,...
                [0:stepSize:nanmax(P(:)); ones(1,[nanmax(P(:))./stepSize]+1).*nanmax(P(:))]);
        end
        if nargin < 4 || isempty(mask)
            mask ={true(size(HD))};
        end

        P = P(:,mask{1});
        HD = HD(:,mask{1});
        T = T(:,mask{1});
        
        wXp = bsxfun(@minus,walls,permute(P,[1 4 3 2]));

        rmat = repmat([permute(cosd(HD),[1 3 4 2]) permute(-sind(HD),[1 3 4 2]); ...
            permute(sind(HD),[1 3 4 2]) permute(cosd(HD),[1 3 4 2])],[1 1 length(walls(1,1,:)) 1]);


        wXphd = nan(size(wXp));
        t = nan([length(walls(1,1,:)) length(walls(1,:,1)) length(HD)]);
        d = nan([length(walls(1,1,:)) length(walls(1,:,1)) length(HD)]);
        for i = 1:length(HD)
            for j = 1:length(walls(1,1,:))
                tmp = rmat(:,:,j,i)*wXp(:,:,j,i);
                [t(j,:,i) d(j,:,i)] = cart2pol(tmp(1,:),tmp(2,:));
            end
        end

        d = floor(d./binDist)+1;
        t = floor(mod(t,2.*pi)./binAng)+1;

        uLoc = unique([d(:) t(:)],'rows');

    %     uLoc(uLoc(:,1)>(nanmax(d(:))./2),:) = [];
        uLoc(uLoc(:,1)>30./binDist,:) = [];
        
        
        prepStack.d = d;
        prepStack.t = t;
        prepStack.uLoc = uLoc;
        prepStack.walls = walls;
    end
    
    maps = nan([length(uLoc(:,2)) length(T(:,1))]);
    smap = nan([length(uLoc(:,2)) 1]);
    parfor i = 1:length(uLoc(:,1))
        isGood = permute(any(any((uLoc(i,1)==d & uLoc(i,2) == t),2),1),[2 3 1 4]);
        smap(i) = nansum(isGood,2);
        maps(i,:) = nansum(T(:,isGood),2); %nansum(bsxfun(@times,isGood,T),2);
    end
    
    s1 = nan(2.*pi./binAng,30./binDist);
    m1 = nan(2.*pi./binAng,30./binDist,length(T(:,1)));
    for i = 1:length(uLoc(:,1))
        s1(uLoc(i,2),uLoc(i,1)) = smap(i);
        m1(uLoc(i,2),uLoc(i,1),:) = maps(i,:);
    end
    smap = s1;
    maps = m1;
    
    tsmap = nansum(smap,4);
%     tsmap(tsmap<300) = nan;
    unm = bsxfun(@rdivide,nansum(maps,4),tsmap);
    m = unm;
    bad = isnan(m);
    m(bad) = 0;    
    m = imfilter(repmat(m,[3 1 1 1]),fspecial('gauss',round([5.*kern 5.*kern]),kern),'same');
    m = m([length(m(:,1,1,1)).*(1./3)]+1:[length(m(:,1,1,1)).*(2./3)],:,:,:);
    m(bad) = nan;
end






























