function [rate] = FiringRateSortAmySocial(shuf,itterations,ms)
% %{

rate = zeros(length(shuf.mrallO),1);

for i = 1 : length(shuf.mrallO)/itterations
    rate(i:length(shuf.mrallO)/itterations:length(shuf.mrallO),1) = sum(ms.deconvolvedSig(:,i))/length(ms.deconvolvedSig(:,1))*30; %length(find(shuf.firing(:,i)))/length(shuf.firing(:,1));
end