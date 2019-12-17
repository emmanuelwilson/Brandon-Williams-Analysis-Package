%sorts shuffling data by mean firing
indOut = [];%sort(indOut9,'descend');
ms = ms;
shuffle = out.mrall;

binarized = Binarize(ms);
btraces = binarized.binarizedTraces;
FiringMean = zeros(length(btraces(1,:)),1);

for i = 1 : length(btraces(1,:))
    FiringMean(i) = length(find(btraces(:,i)))/length(btraces(:,1));
end
Fmean = FiringMean;
for i = 1: length(indOut)   
    Fmean(indOut(i)) = [];  
end

scatmat = NaN(length(shuffle),2);

for i = 1: length(Fmean(:,1))
    scatmat(i:length(Fmean(:,1)):length(shuffle(:,1)),1) = Fmean(i,1);
    scatmat(i:length(Fmean(:,1)):length(shuffle(:,1)),2) = shuffle(i:length(Fmean(:,1)):length(shuffle(:,1)));
end

hold on
figure(1)
scatter(scatmat(:,1),scatmat(:,2))
xlabel('Mean Firing (Active frames/total frames)')
ylabel('Shuffled MRL values')
title('Shuffled MRL vs Mean Firing Rate')