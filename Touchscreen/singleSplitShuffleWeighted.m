%%SplitHalf reliability with random sorting of trial order 100 times.
%
%INPUT:
%   -trial: n X m X p, calcium traces. Where n is the cell number, m is
%   frames, p the number of trials.
%   -cl: Vector containing different trial lengths
%   -cut:Index of what kind of anchoring, 1 for Delay period, 2 for Front
%   Anchored, 3 for Front Anchored. Any excess transient is ignored.
%OUTPUT:
%   -out: Correlation matrix for each split half itteration.

function [out] = singleSplitShuffleWeighted(trial)
if ~isempty(trial)
    if ~(length(trial(1,1,:)) == 0)
        out = zeros(length(trial(:,1,1)),1);
        splitcorr = zeros(length(trial(:,1,1)),100);              
        if ~isempty(trial) && length(trial(1,1,:))>1
            temp = permute(trial,[3 2 1]);
            valsum = sum(temp,2);            
            for i = 1 : length(temp(1,1,:))
                damp = length(find(valsum(:,1,i)))/length(temp(:,1,1));
                if i == 26
                end
                for j = 1: 100
                    rearange = randperm(length(temp(:,1,1)));
                    h1 = nanmean(temp(rearange(1:round(length(rearange)/2)),:,i),1);
                    h2 = nanmean(temp(rearange(round(length(rearange)/2)+1:end),:,i),1);
                    splitcorr(i,j) = damp*corr2(h1,h2);
                    if isnan(splitcorr(i,j))
                        splitcorr(i,j) = 0;
                    end
                end
            end
            out = nanmean(splitcorr,2);
        end
    else
        out = NaN;
    end
else
    out = NaN;
end
end
