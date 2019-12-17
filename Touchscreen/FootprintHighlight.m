%% Footprint Highlight
%Change the wanted cell through "cellNum"
%Modify mask size on line 6 (default 90)
cellNum = 2;
mask = ms.SFPs(:,:,cellNum);
maskThresh = prctile(mask(find(ms.SFPs(:,:,cellNum))),90);
maskind = find(mask>=maskThresh);
mask = zeros(size(mask));
mask(maskind) = 1;
image(imfuse(max(ms.SFPs,[],3),10*mask));