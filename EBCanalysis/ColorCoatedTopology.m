%% Will colour coat the spatial footprints based on preferred firing distance or angle.
% INPUT: 
%       - distanceMat: Matrix of preferred firing distances
%       - angleMat: Matrix of preferred firing angles
%       - ms.SFPs: Miniscope structure, must containt "SFPs" field
%       containing the SpatialFootPrints. 
% OUTPUT:
%       Produces colourcoated figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Émmanuel Wilson

distance = false;                                                          % Change to true if you wish to produce distance based figure
totaldist = 30;
goodcells = ms.SFPs;
if distance
    pref = distanceMat;
else
    pref = angleMat;
end
goodcells(:,:,noPassed) = [];
mask = zeros(size(goodcells(:,:,:)));
for i = 1 : length(goodcells(1,1,:))
    mtemp = goodcells(:,:,i);
    maskThresh = prctile(mtemp(find(mtemp)),90);
    mind = find(mtemp>=maskThresh);    
    mtemp(mind) = pref(i);
    mtemp(find(mtemp~=pref(i))) = 0;
    mask(:,:,i) = mtemp;
end

mask = sum(mask,3);
mask(find(mask == 0)) = NaN;
figure
imagesc(mask,'AlphaData',~isnan(mask))
if distance
    colormap(jet)
    caxis([0 totaldist])
else
    colormap(hsv(360))
    caxis([0 360])
end
colorbar
set(gca,'color',0*[1 1 1]); 
title('5JH2 ContextB Preferred Angle, 2 Changepoints ')