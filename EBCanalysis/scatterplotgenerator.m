% binarized1 = Binarize(ms1);
% binarized2 = Binarize(ms2);
% binarized3 = Binarize(ms3);
% binarized4 = Binarize(ms4);
% 
% [~,~,ind1] = CellFiltering(binarized1.binarizedTraces,0.5);
% [~,~,ind2] = CellFiltering(binarized2.binarizedTraces,0.5);
% [~,~,ind3] = CellFiltering(binarized3.binarizedTraces,0.5);
% [~,~,ind4] = CellFiltering(binarized4.binarizedTraces,0.5);
% ind2 = ind2 + 526;
% ind3 = ind3 + 647;
% ind4 = ind4 + 878;
% indOut = [ind1; ind2; ind3; ind4; dup10j; dupli];
% mrlout = find(mrl<0.1);
% meanout = find(meanF<0.5);
% indOut = cat(1, indOut, mrlout');
% indOut = cat(1, indOut, meanout);
% indOut = unique(indOut);
% moreout1 = [71;72;130;137;139;150;155;158;195;199;217;222;224;264;292;338;347;348;349;416;425;438;443;445;498;503;517;525;526];
% moreout3 =[799;824];
% moreout4 = [899;925;931;933;936;937;938;948;950;959;963;968;972;977;979;980];
% indOut = cat(1,indOut,moreout1);
% indOut = cat(1,indOut,moreout1);
% indOut = cat(1,indOut,moreout3);
% indOut = cat(1,indOut,moreout4);
% indOut = unique(indOut);

% testsnr = snr;
% testMaxPeak = distPeak;
testmrl = out.mrall;
testrm = out.rm;
% testcorr = corr;
% testebc = EBCcells;
testmeanRM = meanF;
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
testrm(:,:,indOut) = [];
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
distanceMat(i,1) = distance;
angle = mean(columns);
angleMat(i,1) = angle;
% pause
clf
end
%
angleup = find(angleMat>360);
angledown = find(angleMat<0);
angleMat(angleup) = angleMat(angleup)-360;
angleMat = angleMat-180;
angledown = find(angleMat<0);
angleMat(angledown) = angleMat(angledown) +360;

figure
histogram(testMaxPeak,round(sqrt(length(testMaxPeak))))
xlabel('Max Peak Value (spike/occ)')
ylabel('Number of Cells')
title('Max Peak Distribution')

figure
histogram(testsnr,round(sqrt(length(testsnr))))
xlabel('Signal to Noise Ratio (Max/Mean)')
ylabel('Number of Cells')
title('SNR Distribution')

figure
histogram(testmrl,round(sqrt(length(testmrl))))
xlabel('Mean Resultant Length')
ylabel('Number of Cells')
title('MRL Distribution')

figure
histogram(testcorr,round(sqrt(length(testcorr))))
xlabel('Correlation Coefficient between halves')
ylabel('Number of Cells')
title('Correlation Distribution')

figure
scatter(testcorr,testsnr)
hold on
scatter(testcorr(testebc),testsnr(testebc),'r','filled')
xlabel('Signal to Noise Ratio (Max/Mean)')
ylabel('Correlation Coefficient between halves')
title('Splithalf Correlation vs SNR')

figure
scatter(testcorr,testMaxPeak)
hold on
scatter(testcorr(testebc),testMaxPeak(testebc),'r','filled')
ylabel('Max Peak Value (spike/occ)')
xlabel('Correlation Coefficient between halves')
title('Splithalf Correlation vs Max Value')

figure
scatter(testcorr,testmrl)
hold on
scatter(testcorr(testebc),testmrl(testebc),'r','filled')
ylabel('Mean Resultant Length')
xlabel('Correlation Coefficient between halves')
title('Splithalf Correlation vs MRL')

figure
scatter(testmrl,testsnr)
hold on
scatter(testmrl(testebc),testsnr(testebc),'r','filled')
ylabel('Signal To Noise Ratio(max/mean)')
xlabel('MRL')
title('MRL vs. SNR')

figure
scatter(testMaxPeak,testsnr)
hold on
scatter(testMaxPeak(testebc),testsnr(testebc),'r','filled')
ylabel('Signal To Noise Ratio(max/mean)')
xlabel('Max Peak Value (spike/occ)')
title('Max Peak vs. SNR')

figure
scatter(testmrl,testMaxPeak)
hold on
scatter(testmrl(testebc),testMaxPeak(testebc),'r','filled')
ylabel('Max Peak Value(spike/occ)')
xlabel('MRL')
title('MRL vs Max Peak Value')

figure
scatter3(testmrl,testMaxPeak, testsnr)
hold on
scatter3(testmrl(testebc),testMaxPeak(testebc), testsnr(testebc),'filled','r')
xlabel('Mean Resultant Length')
ylabel('Max Peak Value (spike/occ)')
zlabel('Signal to noise ratio (max/mean)')
title('3D-Scatter of MLR vs Max Peak vs SNR')

figure
scatter3(testmrl,testMaxPeak, testcorr)
hold on
scatter3(testmrl(testebc),testMaxPeak(testebc), testcorr(testebc),'filled','r')
xlabel('Mean Resultant Length')
ylabel('Max Peak Value (spike/occ)')
zlabel('Splithalf Correlation')
title('3D-Scatter of MLR vs Max Peak vs Correlation')

figure
scatter3(testmrl,testcorr, testsnr)
hold on
scatter3(testmrl(testebc),testcorr(testebc), testsnr(testebc),'filled','r')
xlabel('Mean Resultant Length')
ylabel('Splithalf Correlation')
zlabel('Signal to noise ratio (max/mean)')
title('3D-Scatter of MLR vs Correlation vs SNR')

figure
scatter3(testcorr,testMaxPeak, testsnr)
hold on
scatter3(testcorr(testebc),testMaxPeak(testebc), testsnr(testebc),'filled','r')
xlabel('Splithalf Correlation')
ylabel('Max Peak Value (spike/occ)')
zlabel('Signal to noise ratio (max/mean)')
title('3D-Scatter of Correlation vs Max Peak vs SNR')
