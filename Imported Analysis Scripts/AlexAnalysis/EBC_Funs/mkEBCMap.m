function [m wm unm unwm smap] = mkEBCMap(P,HD,T,masks,walls,prepStack)
    
    
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

        if nargin < 5 || isempty(walls)
            walls = [0 0 0 nanmax(P(:)); 0 0 nanmax(P(:)) 0; ...
                    0 nanmax(P(:)) nanmax(P(:)) nanmax(P(:)); ...
                    nanmax(P(:)) 0 nanmax(P(:)) nanmax(P(:))];
            walls = cat(3,[zeros(1,[nanmax(P(:))./stepSize]+1); 0:stepSize:nanmax(P(:))],...
                [0:stepSize:nanmax(P(:)); zeros(1,[nanmax(P(:))./stepSize]+1)], ...
                [ones(1,[nanmax(P(:))./stepSize]+1).*nanmax(P(:)); 0:stepSize:nanmax(P(:))] ,...
                [0:stepSize:nanmax(P(:)); ones(1,[nanmax(P(:))./stepSize]+1).*nanmax(P(:))]);
        end
        if nargin < 4 || isempty(masks)
            masks ={true(size(HD))};
        end

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
        
        uLoc(uLoc(:,1)>30./binDist,:) = [];
%         uLoc(uLoc(:,1)>14,:) = [];
        
        prepStack.d = d;
        prepStack.t = t;
        prepStack.uLoc = uLoc;
        prepStack.walls = walls;
    end
    
    
%     maps = nan([nanmax(uLoc(:,2)) nanmax(uLoc(:,1)) length(T(:,1)) length(walls(1,1,:))]);
%     smap = nan([nanmax(uLoc(:,2)) nanmax(uLoc(:,1)) 1 length(walls(1,1,:))]);
    
    maps = nan([length(uLoc(:,2)) length(T(:,1)) length(walls(1,1,:))]);
    smap = nan([length(uLoc(:,2)) 1 length(walls(1,1,:))]);

    rT = repmat(T,[1 1 length(walls(1,1,:))]);
    parfor i = 1:length(uLoc(:,1))
        isGood = permute(any((uLoc(i,1)==d & uLoc(i,2) == t),2),[2 3 1 4]);
        smap(i,1,:) = nansum(isGood,2);
        maps(i,:,:) = permute(nansum(bsxfun(@times,isGood,rT),2),[ 2 4 1 3]);
    end
    
    s1 = nan(2.*pi./binAng,30./binDist,length(maps(1,1,:)));
    m1 = nan(2.*pi./binAng,30./binDist,length(T(:,1)),length(maps(1,1,:)));
    for wi = 1:length(maps(1,1,:))
        for i = 1:length(uLoc(:,1))
            s1(uLoc(i,2),uLoc(i,1),wi) = smap(i,1,wi);
            m1(uLoc(i,2),uLoc(i,1),:,wi) = maps(i,:,wi);
        end
    end
    smap = s1;
    maps = m1;
    
    

    tsmap = nansum(smap,3);
%     tsmap(tsmap<300);
    unm = bsxfun(@rdivide,nansum(maps,4),tsmap);
    m = unm;
    bad = isnan(m);
    m(bad) = 0;    
    m = imfilter(repmat(m,[3 1 1 1]),fspecial('gauss',round([5.*kern 5.*kern]),kern),'same');
    m = m([length(m(:,1,1,1)).*(1./3)]+1:[length(m(:,1,1,1)).*(2./3)],:,:,:);
    m(bad) = nan;
    
    for wi = 1:length(walls(1,1,:))
        tsmap = smap(:,:,wi);
%         tsmap(tsmap<120);
        wm = bsxfun(@rdivide,maps(:,:,:,wi),tsmap);
        bad = isnan(wm);
        wm(bad) = 0;
        wm = imfilter(repmat(wm,[3 1 1 1]),fspecial('gauss',round([5.*kern 5.*kern]),kern),'same');
        wm = wm([length(wm(:,1,1,1)).*(1./3)]+1:[length(wm(:,1,1,1)).*(2./3)],:,:,:);
        wm(bad) = nan;
        twm(:,:,:,wi) = wm;
    end
    
    wm = twm;

    
%     tsmap = nansum(smap,3);
%     tsmap = imfilter(repmat(tsmap,[3 1 1 1]),fspecial('gauss',round([5.*kern 5.*kern]),kern),'same');
%     tsmap = tsmap([length(tsmap(:,1,1,1)).*(1./3)]+1:[length(tsmap(:,1,1,1)).*(2./3)],:,:,:);
%     
%     m = nansum(maps,4);
%     bad = isnan(m);
%     m(bad) = 0;
%     m = imfilter(repmat(m,[3 1 1 1]),fspecial('gauss',round([5.*kern 5.*kern]),kern),'same');
%     m = m([length(m(:,1,1,1)).*(1./3)]+1:[length(m(:,1,1,1)).*(2./3)],:,:,:);
%     
%     m = bsxfun(@rdivide,m,tsmap);
%     
%     for wi = 1:length(walls(1,1,:))
%         tsmap = smap(:,:,wi);
%         tsmap = imfilter(repmat(tsmap,[3 1 1 1]),fspecial('gauss',round([5.*kern 5.*kern]),kern),'same');
%         tsmap = tsmap([length(tsmap(:,1,1,1)).*(1./3)]+1:[length(tsmap(:,1,1,1)).*(2./3)],:,:,:);
% 
%         wm = maps(:,:,:,wi);
%         bad = isnan(wm);
%         wm(bad) = 0;
%         wm = imfilter(repmat(wm,[3 1 1 1]),fspecial('gauss',round([5.*kern 5.*kern]),kern),'same');
%         wm = wm([length(wm(:,1,1,1)).*(1./3)]+1:[length(wm(:,1,1,1)).*(2./3)],:,:,:);
% 
%         twm(:,:,:,wi) = bsxfun(@rdivide,wm,tsmap);
%     end
%     
%     wm = twm;

end






























