%% Serialize correlation values from rotated correlation values

function [meanval] = rmcorrserialize(vals)
tempcat = [];
for i = 1 : length(vals(:,1,1))    
    for j = 1 : length(vals(1,:,1))-1
        if any(~isnan(vals(i,j,:)))            
            tempcat = cat(1,tempcat,permute(vals(i,j,~isnan(vals(i,j,:))),[3 2 1]));
        end
    end
end

meanval = tempcat;

end