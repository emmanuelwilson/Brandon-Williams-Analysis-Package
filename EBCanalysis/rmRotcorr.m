%% Rotate Ego RM and correlate, finds correlation dist

function [cordist] = rmRotcorr(ebc1, ebc2,Singlemap,mapind,id1,id2)

c1 = Singlemap(mapind,id1);
c2 = Singlemap(mapind,id2);

cordist = zeros(length(ebc1.rm(1,:,1)),length(c1));

if length(ebc1.rm(:,1,1)) > length(ebc2.rm(:,1,1))
    ebc1.rm = ebc1.rm(1:length(ebc2.rm(:,1,1)),:,:);
elseif length(ebc1.rm(:,1,1)) < length(ebc2.rm(:,1,1))
    ebc2.rm = ebc2.rm(1:length(ebc1.rm(:,1,1)),:,:);
end

for i = 1 : length(c1)
    cordist(1,i) = corr2(ebc1.rm(:,:,c1(i)),ebc2.rm(:,:,c2(i)));
end

for i = 1 : length(cordist(:,1))-1
    rmtemp = circshift(ebc1.rm,i,2);
    for j  = 1 : length(cordist(1,:))
        cordist(i+1,j) = corr2(rmtemp(:,:,c1(j)),ebc2.rm(:,:,c2(j)));
    end
end
end