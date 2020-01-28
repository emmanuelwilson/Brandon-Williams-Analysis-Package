%% average across ratemap correlation values, across all values, outputs single mean value with std error
%takes m x n x n matrix and outputs a single value
function [meanval, valstder] = rmvaltotalmean(vals)
tempcat = [];
for i = 1 : length(vals(:,1,1))    
    for j = 1 : length(vals(1,:,1))-1
        if any(~isnan(vals(i,j,:)))            
            tempcat = cat(1,tempcat,permute(vals(i,j,~isnan(vals(i,j,:))),[3 2 1]));
        end
    end
end

meanval = nanmean(tempcat);
valstder = std(tempcat)/sqrt(length(tempcat));

end