%% takes the normalized mean of the population level calcium
% INPUT: 
%   -minmax: contains the minimum and maximum value of each cell
%   -calcium: NxMxT matrix of calcium transient, N = cell, M = frame/time,
%       T = trial
% OUTPUT:
%    -out: mean of calcium across all cells and trials. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Émmanuel Wilson 

function out = avgfluo(minmax,calcium)
    for i = 1 : length(calcium(:,1,1))
        for j = 1 : length(calcium(1,1,:))
            calcium(i,:,j) = (calcium(i,:,j)-minmax(i,1))./(minmax(i,2) - minmax(i,1));
        end
    end
    calcium = nanmean(calcium,3);
    out = nanmean(calcium,1);
end