%% Rotate Ego RM and correlate, finds correlation dist, single cell

function [cordist] = rmRotcorr_singlecell(rm1, rm2)

if length(rm1(:,1,1)) > length(rm2(:,1,1))
    rm1 = rm1(1:length(rm2(:,1,1)),:,:);
elseif length(rm1(:,1,1)) < length(rm2(:,1,1))
    rm2 = rm2(1:length(rm1(:,1,1)),:,:);
end
cordist = zeros(length(rm1(1,:)),1);
cordist(1) = corr2(rm1,rm2);

for i = 1 : length(cordist(:,1))-1
    rmtemp = circshift(rm1,i,2);
    cordist(i+1,1) = corr2(rmtemp(:,:),rm2(:,:));
end
end