%% average across ratemap correlation values 
%takes m x n x n matrix and outputs an m x 1 matrix
function [meanval, valstder] = rmvalmean(vals)

for i = 1 : length(vals(:,1,1))
    tempcat = [];
    for j = 1 : length(vals(1,:,1))-1
        if any(~isnan(vals(i,j,:)))            
            tempcat = cat(1,tempcat,permute(vals(i,j,~isnan(vals(i,j,:))),[3 2 1]));
        end
    end
    meanval(i,1) = nanmean(tempcat);
    valstder(i,1) = std(tempcat)/sqrt(length(tempcat));
end

end