function [cellmap,exclude] = msExtractSFPsCellReg2020(SFP)
%msExtractSFPs Extracts spatial footprints to perform chronic re-alignment
% Converts spatial footprints from m,k,n (UCLA) to n,m,k (Ziv's lab) where
% n is the number of neurons, k pixels in x axis, and m pixels in y axis
%
% Author: Guillaume Etter
% Contact: etterguillaume@gmail.com
% modified: Emmanuel Wilson

exclude = ones(1,size(SFP,3));
SFP_temp = permute(sum(sum(SFP,1),2),[3 2 1]);
exclude(find(SFP_temp == 0)) = 0;
exclude = logical(exclude);
cellmap = find(SFP_temp);

end