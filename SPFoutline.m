function [ms] = SPFoutline(ms)
BW = zeros(size(ms.SFPs));
for cellNum = 1: length(ms.SFPs(1,1,:))
    mask = ms.SFPs(:,:,cellNum);
    maskThresh = prctile(mask(find(ms.SFPs(:,:,cellNum))),90);
    maskind = find(mask>=maskThresh);
    mask = zeros(size(mask));
    mask(maskind) = 1;
    outline = bwboundaries(mask);
    if length(outline) == 1
        for i = 1 : length(outline{1}(:,1))
            BW(outline{1}(i,1),outline{1}(i,2),cellNum) = 1;
        end
    else
       fprintf(['Cell ', num2str(cellNum), ' is identified as ', num2str(length(outline)),' objects',newline]);
    end
    outline = [];
end

ms.outlines = sum(BW,3);
end