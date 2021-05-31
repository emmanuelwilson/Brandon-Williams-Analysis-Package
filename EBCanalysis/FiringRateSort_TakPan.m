function [rate] = FiringRateSort_TakPan(shuf,itterations)

rate = zeros(length(shuf.mrall),1);

for i = 1 : length(shuf.mrall)/itterations
    rate(i:length(shuf.mrall)/itterations:length(shuf.mrall),1) = sum(shuf.firing(:,i))/length(shuf.firing(:,1))*5;
end

end