testrm = out.rm;
% testcorr = corr;
% testebc = EBCcells;
% testmeanRM = meanF;
cellout = [];

distanceMat = [];
angleMat = [];

% for i =1 : length(testebc)
%     if ~isempty(find(testebc(i)== indOut))
%         cellout = i;
%     else        
%         testebc(i) = testebc(i) - sum((indOut<testebc(i)));
%     end
% end

if ~isempty(cellout)
%     testebc(cellout) = [];
end

% testsnr(indOut) = [];
% testMaxPeak(indOut) = [];
% testmrl(indOut) = [];
% testmeanRM(indOut) = [];
testrm(:,:,noPassed7Ap) = [];
% testcorr(indOut) = [];

%Mark's code
nCells = size(testrm,3);
BWtot = zeros(length(testrm(:,1,1)),722);
for i = 1:nCells
rmcrit = testrm(:,:,i);
rmcrit = [rmcrit rmcrit];
% if max(max(rmcrit))<1, continue ,end
ind = find(rmcrit<max(max(rmcrit))*.5);
rm_thresh = rmcrit;
rm_thresh(ind) = 0;
BW = imregionalmax(rm_thresh,26);
% figure(1), 
% subplot(3,1,1), imagesc(rmcrit)
% xlabel('Angle (degrees)')
% ylabel('Distance (cm)')
% subplot(3,1,2), imagesc(rm_thresh)
% xlabel('Angle (degrees)')
% ylabel('Distance (cm)')
BW2 = bwareafilt(rm_thresh>0,1,'largest');
BWtot = BWtot + BW2;
imagesc(BWtot)
xlabel('Angle (degrees)')
ylabel('Distance (cm)')
[rows, columns] = find(BW2>0);
distance = mean(rows);
distanceMat(i,1) = distance;
angle = mean(columns);
angleMat(i,1) = angle;
pause(0.01)
end
%
angleup = find(angleMat>360);
angledown = find(angleMat<0);
angleMat(angleup) = angleMat(angleup)-360;
angleMat = angleMat-180;
angledown = find(angleMat<0);
angleMat(angledown) = angleMat(angledown) +360;