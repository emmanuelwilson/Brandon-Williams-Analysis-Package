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

distance = 1;                                                          % Change to true if you wish to produce distance based figure
totaldist = 35;
goodcells = ms104.SFPs(:,:,fullpass104);
mouse = '10JH3';
con = 'A4';
if distance
    pref = distanceMat104;
    sname = 'Distance';
else
    pref = angleMat104;
    sname = 'Angle';
end
% % % goodcells(:,:,noPassed) = [];
mask = zeros(size(goodcells(:,:,:)));
for i = 1 : length(goodcells(1,1,:))
    mtemp = goodcells(:,:,i);
    maskThresh = prctile(mtemp(find(mtemp)),90);
    mind = find(mtemp>=maskThresh);    
    mtemp(mind) = pref(i);
    mtemp(find(mtemp~=pref(i))) = 0;
    mask(:,:,i) = mtemp;
%     BW(:,:,i) = bwperim(mtemp);
end

% BW = max(BW,[],3);
mask = max(mask,[],3);
mask(find(mask == 0)) = NaN;
% mask(find(BW)) = 365;
figure(1)
imagesc(mask,'AlphaData',~isnan(mask))
if distance
    t = jet(255);
%     t = cat(1,t,[0 0 0]);
    colormap(t)
    caxis([0 totaldist])
else    
    t = hsv(255);
%     t = cat(1,t,[0 0 0]);    
    colormap(t)
    caxis([0 361])    
end
colorbar
set(gca,'color',[0 0 0]); 
title([mouse 'Context' con ' Preferred Firing ' sname],'FontSize',26)
% 
saveas(gcf,['Pref', sname,'TopologyContext', con,'_', mouse, '.fig'])
saveas(gcf,['Pref', sname,'TopologyContext ', con, '_', mouse, '.eps'])