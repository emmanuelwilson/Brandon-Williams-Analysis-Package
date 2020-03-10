%% Recreates EBC ratemaps and identifies the angle and distance of the largest object.
% INPUT: 
%       -out.rm: the EBC output which contains the rate map variable "rm"
%       which is a n X m X p matrix. Where n and m are the angle and distance
%       dimensions of the ratemap and p is the cell#
%       -noPassed: Cell index of cells you do not wish to include in the
%       analysis.
% OUTPUT:
%       -distanceMat: Matrix of preffered firing distance
%       -angleMat: Matrix of preffered firing angle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Mark Brandon/Émmanuel Wilson

anglebin = 6;
distancebin = 1.2;

testrm = EBCA14.rm;
% testcorr = corr;
% testebc = EBCcells;
% testmeanRM = meanF;
cellout = [];

distanceMatA14 = [];
angleMatA14 = [];

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
% testrm(:,:,noPassed) = [];%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% testcorr(indOut) = [];

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
subplot(3,1,1), imagesc(rmcrit)
% title(num2str(passed9(i)))%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlabel('Angle (degrees)')
ylabel('Distance (cm)')
subplot(3,1,2), imagesc(rm_thresh)
xlabel('Angle (degrees)')
ylabel('Distance (cm)')
BW2 = bwareafilt(rm_thresh>0,1,'largest');
subplot(3,1,3), imagesc(BW2)
xlabel('Angle (degrees)')
ylabel('Distance (cm)')
[rows, columns] = find(BW2>0);
distance = mean(rows);
distanceMatA14(i,1) = distance;
angle = mean(columns);
angleMatA14(i,1) = angle;
vline(angleMatA14(i,1),'w')
hline(distanceMatA14(i,1),'w')
% pause
clf
end
%
distanceMatA14 = distanceMatA14*distancebin;
angleMatA14 = angleMatA14 *anglebin;
angleup = find(angleMatA14>360);
angledown = find(angleMatA14<0);
angleMatA14(angleup) = angleMatA14(angleup)-360;
angleMatA14 = angleMatA14-180;
angledown = find(angleMatA14<0);
angleMatA14(angledown) = angleMatA14(angledown) +360;