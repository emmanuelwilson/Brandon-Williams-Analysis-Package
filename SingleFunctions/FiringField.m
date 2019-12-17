testrm = out.rm;

distanceMat = [];
angleMat = [];
Ffield = zeros(length(out.rm(:,1,1)),722);

testrm(:,:,indOut) = [];


%Mark's code
nCells = size(testrm,3);
for i = 1:nCells
rmcrit = testrm(:,:,i);
rmcrit = [rmcrit rmcrit];
% if max(max(rmcrit))<1, continue ,end
ind = find(rmcrit<max(max(rmcrit))*.5);
rm_thresh = rmcrit;
rm_thresh(ind) = 0;
BW = imregionalmax(rm_thresh,26);
figure(1), 
subplot(4,1,1), imagesc(rmcrit)
xlabel('Angle (degrees)')
ylabel('Distance (cm)')
subplot(4,1,2), imagesc(rm_thresh)
xlabel('Angle (degrees)')
ylabel('Distance (cm)')
BW2 = bwareafilt(rm_thresh>0,1,'largest');
subplot(4,1,3), imagesc(BW2)
Ffield = Ffield + BW2;
subplot(4,1,4), imagesc(Ffield)
xlabel('Angle (degrees)')
ylabel('Distance (cm)')
[rows, columns] = find(BW2>0);
% pause
% clf
end

Fbig = Ffield(:,360:end-3);
Ffield(:,361:end) = [];
Ffield = Ffield + Fbig;

J = find(Ffield == 0);
imAlpha=ones(size(Ffield));
imAlpha(J)=0;

figure
imagesc(Ffield,'AlphaData',imAlpha);
set(gca,'color',[1 1 1]);
xlabel('Angle (deg)')
ylabel('Distance (cm)')
title('Total Firing Field')
colorbar
