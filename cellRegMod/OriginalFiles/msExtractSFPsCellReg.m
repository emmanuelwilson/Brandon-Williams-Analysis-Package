function [SFP,cellmap] = msExtractSFPsCellReg(ms)
%msExtractSFPs Extracts spatial footprints to perform chronic re-alignment
% Converts spatial footprints from m,k,n (UCLA) to n,m,k (Ziv's lab) where
% n is the number of neurons, k pixels in x axis, and m pixels in y axis
%
% Author: Guillaume Etter
% Contact: etterguillaume@gmail.com
% modified: Emmanuel Wilson

skipcount = 0;
cellmap = 1: length(ms.SFPs(1,1,:));
for cell_i = 1:size(ms.SFPs,3)
    SFP_temp = ms.SFPs(:,:,cell_i);
    SFP_temp(SFP_temp<0.5*max(max(SFP_temp))) = 0; % This is to sharpen footprints, based on Ziv lab method
    if sum(sum(SFP_temp)) > 0
        SFP(cell_i-skipcount,:,:) = SFP_temp;
    else
        cellmap(find(cell_i)) = [];
        skipcount = skipcount +1;
    end
end

end

