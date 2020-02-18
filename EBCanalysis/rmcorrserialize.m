%% Serialize correlation values from rotated correlation values

function [meanval] = rmcorrserialize(vals)
tempcat = [];
t = 0;
for i = 1 : length(vals(:,1,1))    
    for j = 1 : length(vals(1,:,1))
        if any(~isnan(vals(i,j,:)))
            tempcat = cat(1,tempcat,permute(vals(i,j,~isnan(vals(i,j,:))),[3 2 1]));
            t = 1;
        end
    end
    if ~t
        testytest =  permute(vals(i,:,:),[3 2 1]);
    end
    t = 0;
end

meanval = tempcat;

end