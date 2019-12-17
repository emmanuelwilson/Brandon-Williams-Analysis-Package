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

pref1 = 6; % color cel 1

pref2 = 17; % color cel 10

pref3 = 27; % color cel 15

pref4 = 39; % color cel 20

goodcells1 = ms1.SFPs(:,:,g1);
goodcells2 = ms2.SFPs(:,:,g2);
goodcells3 = ms3.SFPs(:,:,g3);
goodcells4 = ms4.SFPs(:,:,g4);

%day 1
mask = zeros(size(goodcells1(:,:,:)));
for i = 1 : length(goodcells1(1,1,:))
    mtemp = goodcells1(:,:,i);
    maskThresh = prctile(mtemp(find(mtemp)),90);
    mind = find(mtemp>=maskThresh);    
    mtemp(mind) = pref1;
    mtemp(find(mtemp~=pref1)) = 0;
    mask(:,:,i) = mtemp;
end

mask = sum(mask,3);
mask1 = mask;
% mask1(find(~mask)) = NaN;

%day 2
mask = zeros(size(goodcells2(:,:,:)));
for i = 1 : length(goodcells2(1,1,:))
    mtemp = goodcells2(:,:,i);
    maskThresh = prctile(mtemp(find(mtemp)),90);
    mind = find(mtemp>=maskThresh);    
    mtemp(mind) = pref2;
    mtemp(find(mtemp~=pref2)) = 0;
    mask(:,:,i) = mtemp;
end

mask = sum(mask,3);
mask2 = mask;
% mask2(find(mask == 0)) = NaN;

%day 3
mask = zeros(size(goodcells3(:,:,:)));
for i = 1 : length(goodcells3(1,1,:))
    mtemp = goodcells3(:,:,i);
    maskThresh = prctile(mtemp(find(mtemp)),90);
    mind = find(mtemp>=maskThresh);    
    mtemp(mind) = pref3;
    mtemp(find(mtemp~=pref3)) = 0;
    mask(:,:,i) = mtemp;
end

mask = sum(mask,3);
mask3 = mask;
% mask3(find(mask == 0)) = NaN;

%day 4
mask = zeros(size(goodcells4(:,:,:)));
for i = 1 : length(goodcells4(1,1,:))
    mtemp = goodcells4(:,:,i);
    maskThresh = prctile(mtemp(find(mtemp)),90);
    mind = find(mtemp>=maskThresh);    
    mtemp(mind) = pref4;
    mtemp(find(mtemp~=pref4)) = 0;
    mask(:,:,i) = mtemp;
end

mask = sum(mask,3);
mask4 = mask;
% mask4(find(mask == 0)) = NaN;

mask = cat(3,mask1,mask2,mask3,mask4);
mask = nansum(mask,3);
mask(find(mask == 0)) = NaN;

figure
imagesc(mask,'AlphaData',~isnan(mask))
caxis([1 89])
colormap(jet)
% caxis([0 totaldist])
colorbar
set(gca,'color',0*[1 1 1]); % color de fondo  0*
title('3101F, 6 seconds Delay, Choice  1')