load('ms.mat');
load('behav.mat');
load('msTouchSync.mat');



%% -------------STEP 3) Reliability ---------------

%Step1: Align data to Anchor points
%Generate raster plots for all trials and all cells
[Calcium, TrialInds] = TouchRaster(events,ms,1,1,1,0,1,0,0);
load('minmax.mat')
load('SepDelay.mat')
load('SepFront.mat')
load('SepBack.mat')

%Step2: Find Reliable Cell
%Cells that pass criteria
ShuffledCrit = TouchScreenShuffle(ms,events,TrialInds);
ShuffledCritSep = TouchScreenShuffleSeperation(SepDelay,SepFront,SepBack,ms,TrialInds,events);
save('ShuffledCrit.mat','ShuffledCrit')

%Figure1: Seqence of population data
%Delay
[Dcorrectmean,rastsort] = sortpeaks(nanmean(ShuffledCrit.DcorrectPassed,3));
Dincorrectmean = nanmean(Calcium.Dincorrect(ShuffledCrit.DcorrectCells,:,:),3);
Dincorrectmean = Dincorrectmean(rastsort,:,:);
Dcorcorrectmean =nanmean(Calcium.DcorrectCor(ShuffledCrit.DcorrectCells,:,:),3);
Dcorcorrectmean = Dcorcorrectmean(rastsort,:,:);
Dcorincorrectmean =nanmean(Calcium.DincorrectCor(ShuffledCrit.DcorrectCells,:,:),3);
Dcorincorrectmean = Dcorincorrectmean(rastsort,:,:);
% Dincorrectmean = sortpeaks(nanmean(Calcium.Dincorrect(ShuffledCrit.DcorrectCells,:,:),3));
% Dcorcorrectmean = sortpeaks(nanmean(Calcium.DcorrectCor(ShuffledCrit.DcorrectCells,:,:),3));
% Dcorincorrectmean = sortpeaks(nanmean(Calcium.DincorrectCor(ShuffledCrit.DcorrectCells,:,:),3));

figure(1)
colormap parula
subplot(2,2,1)
imagesc(Dcorrectmean)
ylabel('Cell')
xlabel('Time (frames)')
title('Delay: Correct Mean Fluorescence')
colorbar
subplot(2,2,2)
imagesc(Dincorrectmean)
ylabel('Cell')
xlabel('Time (frames)')
title('Delay: Incorrect Mean Fluorescence')
colorbar
subplot(2,2,3)
imagesc(Dcorcorrectmean)
ylabel('Cell')
xlabel('Time (frames)')
title('Delay: Correct Correction Mean Fluorescence')
colorbar
subplot(2,2,4)
imagesc(Dcorincorrectmean)
ylabel('Cell')
xlabel('Time (frames)')
title('Delay: Incorrect Correction Mean Fluorescence')
colorbar
saveas(figure(1),'Delay/Delay_PassedCriteria.fig')

if isempty(Dcorincorrectmean) || (nansum(nansum(Dcorincorrectmean))) == 0
    IncCor = false;
else
    IncCor = true;
end
%Front Anchored

[Fcorrectmean,rastsort] = sortpeaks(nanmean(ShuffledCrit.FcorrectPassed,3));
Fincorrectmean = nanmean(Calcium.Fincorrect(ShuffledCrit.FcorrectCells,:,:),3);
Fincorrectmean = Fincorrectmean(rastsort,:,:);
Fcorcorrectmean =nanmean(Calcium.FcorrectCor(ShuffledCrit.FcorrectCells,:,:),3);
Fcorcorrectmean = Fcorcorrectmean(rastsort,:,:);
if ~isempty(Calcium.FincorrectCor)
Fcorincorrectmean =nanmean(Calcium.FincorrectCor(ShuffledCrit.FcorrectCells,:,:),3);
Fcorincorrectmean = Fcorincorrectmean(rastsort,:,:);
else
    Fcorincorrectmean = [];
end
% % Fcorrectmean = sortpeaks(nanmean(ShuffledCrit.FcorrectPassed,3));
% % Fincorrectmean = sortpeaks(nanmean(Calcium.Fincorrect(ShuffledCrit.FcorrectCells,:,:),3));
% % Fcorcorrectmean = sortpeaks(nanmean(Calcium.FcorrectCor(ShuffledCrit.FcorrectCells,:,:),3));
% % if ~isempty(Calcium.FincorrectCor)
% %     Fcorincorrectmean = sortpeaks(nanmean(Calcium.FincorrectCor(ShuffledCrit.FcorrectCells,:,:),3));
% % else
% %     Fcorincorrectmean = [];
% % end

figure(1)
colormap parula
subplot(2,2,1)
imagesc(Fcorrectmean)
ylabel('Cell')
xlabel('Time (frames)')
title('Front Anchored: Correct Mean Fluorescence')
colorbar
subplot(2,2,2)
imagesc(Fincorrectmean)
ylabel('Cell')
xlabel('Time (frames)')
title('Front Anchored: Incorrect Mean Fluorescence')
colorbar
subplot(2,2,3)
imagesc(Fcorcorrectmean)
ylabel('Cell')
xlabel('Time (frames)')
title('Front Anchored: Correct Correction Mean Fluorescence')
colorbar
subplot(2,2,4)
imagesc(Fcorincorrectmean)
ylabel('Cell')
xlabel('Time (frames)')
title('Front Anchored: Incorrect Correction Mean Fluorescence')
colorbar
saveas(figure(1),'FrontAnchored/Front_PassedCriteria.fig')

%Back Anchored
%sort cells by firing order for the correct trials and keep same cell order
%across all trial types
[Bcorrectmean,rastsort] = sortpeaks(nanmean(ShuffledCrit.BcorrectPassed,3));
Bincorrectmean = nanmean(Calcium.Bincorrect(ShuffledCrit.BcorrectCells,:,:),3);
Bincorrectmean = Bincorrectmean(rastsort,:,:);
Bcorcorrectmean =nanmean(Calcium.BcorrectCor(ShuffledCrit.BcorrectCells,:,:),3);
Bcorcorrectmean = Bcorcorrectmean(rastsort,:,:);
if ~isempty(Calcium.BincorrectCor)
Bcorincorrectmean =nanmean(Calcium.BincorrectCor(ShuffledCrit.BcorrectCells,:,:),3);
Bcorincorrectmean = Bcorincorrectmean(rastsort,:,:);
else    
    Bcorincorrectmean = [];   
end
% % Bcorrectmean = sortpeaks(nanmean(ShuffledCrit.BcorrectPassed,3));
% % Bincorrectmean = sortpeaks(nanmean(Calcium.Bincorrect(ShuffledCrit.BcorrectCells,:,:),3));
% % Bcorcorrectmean = sortpeaks(nanmean(Calcium.BcorrectCor(ShuffledCrit.BcorrectCells,:,:),3));
% % if ~isempty(Calcium.BincorrectCor)
% %     Bcorincorrectmean = sortpeaks(nanmean(Calcium.BincorrectCor(ShuffledCrit.BcorrectCells,:,:),3));
% % else
% %     Bcorincorrectmean = [];   
% % end

figure(1)
colormap parula
subplot(2,2,1)
imagesc(Bcorrectmean)
vline(length(Bcorrectmean(1,:)) - 450,'g')
ylabel('Cell')
xlabel('Time (frames)')
title('Back Anchored: Correct Mean Fluorescence')
colorbar
subplot(2,2,2)
imagesc(Bincorrectmean)
vline(length(Bincorrectmean(1,:)) - 450,'r')
ylabel('Cell')
xlabel('Time (frames)')
title('Back Anchored: Incorrect Mean Fluorescence')
colorbar
subplot(2,2,3)
imagesc(Bcorcorrectmean)
vline(length(Bcorcorrectmean(1,:)) - 450,'g')
ylabel('Cell')
xlabel('Time (frames)')
title('Back Anchored: Correct Correction Mean Fluorescence')
colorbar
subplot(2,2,4)
imagesc(Bcorincorrectmean)
ylabel('Cell')
xlabel('Time (frames)')
title('Back Anchored: Incorrect Correction Mean Fluorescence')
colorbar
saveas(figure(1),'BackAnchored/Back_PassedCriteria.fig')

%Figure2: Population Representation of task contingencies
figure(1)
%Delay
Dcorrectpop = avgfluo(minmax(ShuffledCrit.DcorrectCells,:),ShuffledCrit.DcorrectPassed);
Dincorrectpop = avgfluo(minmax(ShuffledCrit.DcorrectCells,:),Calcium.Dincorrect(ShuffledCrit.DcorrectCells,:,:));
Dcorcorrectpop = avgfluo(minmax(ShuffledCrit.DcorrectCells,:),Calcium.DcorrectCor(ShuffledCrit.DcorrectCells,:,:));
if IncCor
Dcorincorrectpop = avgfluo(minmax(ShuffledCrit.DcorrectCells,:),Calcium.DincorrectCor(Shuffled.DcorrectCells,:,:));
else
    Dcorincorrectpop = 0;
end
figure(1)
colormap parula
subplot(2,2,1)
plot(Dcorrectpop)
ylabel('Mean Fluerescence Activity (0-1)')
xlabel('Time (frames)')
title('Delay: Correct Population Level Fluorescence')
subplot(2,2,2)
plot(Dincorrectpop)
ylabel('Mean Fluerescence Activity (0-1)')
xlabel('Time (frames)')
title('Delay: Incorrect Population Level Fluorescence')
subplot(2,2,3)
plot(Dcorcorrectpop)
ylabel('Mean Fluerescence Activity (0-1)')
xlabel('Time (frames)')
title('Delay: Correct Correction Population Level Fluorescence')
subplot(2,2,4)
plot(Dcorincorrectpop)
ylabel('Mean Fluerescence Activity (0-1)')
xlabel('Time (frames)')
title('Delay: Incorrect Correction Population Level Fluorescence')
saveas(figure(1),'Delay/Delay_PassedCriteriaPOPULATION.fig')
%Front Anchored
Fcorrectpop = avgfluo(minmax(ShuffledCrit.FcorrectCells,:),ShuffledCrit.FcorrectPassed);
Fincorrectpop = avgfluo(minmax(ShuffledCrit.FcorrectCells,:),Calcium.Fincorrect(ShuffledCrit.FcorrectCells,:,:));
Fcorcorrectpop = avgfluo(minmax(ShuffledCrit.FcorrectCells,:),Calcium.FcorrectCor(ShuffledCrit.FcorrectCells,:,:));
if IncCor
Fcorincorrectpop = avgfluo(minmax(Shuffled.FcorrectCells,:),Calcium.Fincorrect(Shuffled.FcorrectCells,:,:));
else
    Fcorincorrectpop = 0;
end

figure(1)
colormap parula
subplot(2,2,1)
plot(Fcorrectpop)
ylabel('Mean Fluerescence Activity (0-1)')
xlabel('Time (frames)')
title('Front Anchored: Correct Population Level Fluorescence')
subplot(2,2,2)
plot(Fincorrectpop)
ylabel('Mean Fluerescence Activity (0-1)')
xlabel('Time (frames)')
title('Front Anchored: Incorrect Population Level Fluorescence')
subplot(2,2,3)
plot(Fcorcorrectpop)
ylabel('Mean Fluerescence Activity (0-1)')
xlabel('Time (frames)')
title('Front Anchored: Correct Correction Population Level Fluorescence')
subplot(2,2,4)
plot(Fcorincorrectpop)
ylabel('Mean Fluerescence Activity (0-1)')
xlabel('Time (frames)')
title('Front Anchored: Incorrect Correction Population Level Fluorescence')
saveas(figure(1),'FrontAnchored/Front_PassedCriteriaPOPULATION.fig')

% Back Anchored
Bcorrectpop = avgfluo(minmax(ShuffledCrit.BcorrectCells,:),ShuffledCrit.BcorrectPassed);
Bincorrectpop = avgfluo(minmax(ShuffledCrit.BcorrectCells,:),Calcium.Bincorrect(ShuffledCrit.BcorrectCells,:,:));
Bcorcorrectpop = avgfluo(minmax(ShuffledCrit.BcorrectCells,:),Calcium.BcorrectCor(ShuffledCrit.BcorrectCells,:,:));
if IncCor
    Bcorincorrectpop = avgfluo(minmax(ShuffledCrit.BcorrectCells,:),Calcium.Bincorrect(Shuffled.BcorrectCells,:,:));
else
    Bcorincorrectpop = 0;
end

figure(1)
colormap parula
subplot(2,2,1)
plot(Bcorrectpop)
vline(length(Bcorrectpop) - 450,'g')
ylabel('Mean Fluerescence Activity (0-1)')
xlabel('Time (frames)')
title('Back Anchored: Correct Population Level Fluorescence')
subplot(2,2,2)
plot(Bincorrectpop)
ylabel('Mean Fluerescence Activity (0-1)')
xlabel('Time (frames)')
title('Back Anchored: Incorrect Population Level Fluorescence')
subplot(2,2,3)
plot(Bcorcorrectpop)
vline(length(Bcorcorrectpop) - 450,'g')
ylabel('Mean Fluerescence Activity (0-1)')
xlabel('Time (frames)')
title('Back Anchored: Correct Correction Population Level Fluorescence')
subplot(2,2,4)
plot(Bcorincorrectpop)
ylabel('Mean Fluerescence Activity (0-1)')
xlabel('Time (frames)')
title('Back Anchored: Incorrect Correction Population Level Fluorescence')
saveas(figure(1),'BackAnchored/Back_PassedCriteriaPOPULATION.fig')
clf

%% -----------STEP 4): Memeory related firing ---------------

%Step1: Which cell passed criteria based on choice location
%Figure1: Sample selected activity of individual cells
n = 5;
m = 4;

%Delay
%find all choices
choices = PokeFind(SepDelay);
count = 1;
for i = 1 : length(choices)
    %Correct
    if ~isempty(SepDelay.CorrectTrial{choices(i)})
        passedcorrectD = SepDelay.CorrectTrial{choices(i)}(:,:,ShuffledCrit.DcorrectCells);
        for j = 1 : ShuffledCrit.DcorrectNumberPassed
            cellNum = ShuffledCrit.DcorrectCells(j);
            figure(1)
            subplot(n,m,count)
            count = count +1;
            h = imagesc(passedcorrectD(:,:,j));
            axis('tight')
            colormap('parula')
            colorbar
            caxis([minmax(ShuffledCrit.DcorrectCells(j),1) minmax(ShuffledCrit.DcorrectCells(j),2)])
            view(2)
            xlabel('Time(frame)')
            ylabel('Trial Number')
            title(['Passed Criteria Correct: Cell ',num2str(cellNum)])
            if mod(j,n*m) == 0 || j == ShuffledCrit.DcorrectNumberPassed
                pause(0.01)
                count = 1;                
                saveas(gcf,[pwd,'/Delay/ChoiceSeperation/CorrectTrials/',num2str(choices(i)),'/', num2str(cellNum),'CorrectTrials.jpg']);
                clf
            end
        end
    else
        passedcorrectD = [];
    end
    %Incorrect
    if ~isempty(SepDelay.IncorrectTrial{choices(i)})
        passedincorrectD = SepDelay.IncorrectTrial{choices(i)}(:,:,ShuffledCrit.DincorrectCells);
        for j = 1 : ShuffledCrit.DincorrectNumberPassed
            cellNum = ShuffledCrit.DincorrectCells(j);
            figure(1)
            subplot(n,m,count)
            count = count +1;
            h = imagesc(passedincorrectD(:,:,j));
            axis('tight')
            colormap('parula')
            colorbar
            caxis([minmax(ShuffledCrit.DincorrectCells(j),1) minmax(ShuffledCrit.DincorrectCells(j),2)])
            view(2)
            xlabel('Time(frame)')
            ylabel('Trial Number')
            title(['Passed Criteria Incorrect: Cell ',num2str(cellNum)])
            if mod(j,n*m) == 0 || j == ShuffledCrit.DincorrectNumberPassed
                pause(0.01)
                count = 1;
                saveas(gcf,[pwd,'/Delay/ChoiceSeperation/IncorrectTrials/',num2str(choices(i)),'/', num2str(cellNum),'IncorrectTrials.jpg']);
                clf
            end
        end
    else
        passedincorrectD = [];
    end
    %Correct Correction
    if ~isempty(SepDelay.CorrectCorrectionTrial{choices(i)})
        passedcorrectcorrectionD = SepDelay.CorrectCorrectionTrial{choices(i)}(:,:,ShuffledCrit.DccorCells);
        for j = 1 : ShuffledCrit.DccorNumberPassed
            cellNum = ShuffledCrit.DccorCells(j);
            figure(1)
            subplot(n,m,count)
            count = count +1;
            h = imagesc(passedcorrectcorrectionD(:,:,j));
            axis('tight')
            colormap('parula')
            colorbar
            caxis([minmax(ShuffledCrit.DccorCells(j),1) minmax(ShuffledCrit.DccorCells(j),2)])
            view(2)
            xlabel('Time(frame)')
            ylabel('Trial Number')
            title(['Passed Criteria Correct Correction: Cell ',num2str(cellNum)])
            if mod(j,n*m) == 0 || j == ShuffledCrit.DccorNumberPassed
                pause(0.01)
                count = 1;
                saveas(gcf,[pwd,'/Delay/ChoiceSeperation/CorrectCorrectionTrials/',num2str(choices(i)),'/', num2str(cellNum),'CorrectCorrectionTrials.jpg']);
                clf
            end
        end
    else
        passedcorrectcorrectionD = [];
    end
    %Incorrect Correction
    if IncCor && ~isempty(SepDelay.IncorrectCorrectionTrial{choices(i)})
        passedincorrectcorrectionD = SepDelay.IncorrectCorrectionTrial{choices(i)}(:,:,ShuffledCrit.DicorCells);
        for j = 1 : ShuffledCrit.DicorNumberPassed
            cellNum = ShuffledCrit.DicorCells(j);
            figure(1)
            subplot(n,m,count)
            count = count +1;
            h = imagesc(passedincorrectcorrectionD(:,:,j));
            axis('tight')
            colormap('parula')
            colorbar
            caxis([minmax(ShuffledCrit.DicorCells(j),1) minmax(ShuffledCrit.DicorCells(j),2)])
            view(2)
            xlabel('Time(frame)')
            ylabel('Trial Number')
            title(['Passed Criteria Incorrect Correction: Cell ',num2str(cellNum)])
            if mod(j,n*m) == 0 || j == ShuffledCrit.DicorNumberPassed
                pause(0.01)
                count = 1;
                saveas(gcf,[pwd,'/Delay/ChoiceSeperation/IncorrectCorrectionTrials/',num2str(choices(i)),'/', num2str(cellNum),'IncorrectCorrectionTrials.jpg']);
                clf
            end
        end
    else
        passedincorrectcorrectionD =[];
    end
    meancDpop = sortpeaks(permute(mean(passedcorrectD,1),[3 2 1]));
    meaniDpop = sortpeaks(permute(mean(passedincorrectD,1),[3 2 1]));
    meanccorDpop = sortpeaks(permute(mean(passedcorrectcorrectionD,1),[3 2 1]));
    if IncCor
    meanicorDpop = sortpeaks(permute(mean(passedincorrectcorrectionD,1),[3 2 1]));
    else
        meanicorDpop = [];
    end
    figure(1)
    subplot(2,2,1)
    imagesc(meancDpop)
    ylabel('Cell')
    xlabel('Time (frames)')
    title('Delay: Correct Mean Fluorescence')
    subplot(2,2,2)
    imagesc(meaniDpop)
    ylabel('Cell')
    xlabel('Time (frames)')
    title('Delay: Incorrect Mean Fluorescence')
    subplot(2,2,3)
    imagesc(meanccorDpop)
    ylabel('Cell')
    xlabel('Time (frames)')
    title('Delay: Correct Correction Mean Fluorescence')
    subplot(2,2,4)
    imagesc(meanicorDpop)
    ylabel('Cell')
    xlabel('Time (frames)')
    title('Delay: Incorrect Correction Mean Fluorescence')
    saveas(gcf,[pwd,'/Delay/ChoiceSeperation_',num2str(choices(i)),'_PopulationRasterPlot.jpg']);
    clf    
end

%Front Anchored
choices = PokeFind(SepFront);
count = 1;
for i = 1 : length(choices)
    %Correct
    if ~isempty(SepFront.CorrectTrial{choices(i)})
        passedcorrectF = SepFront.CorrectTrial{choices(i)}(:,:,ShuffledCrit.FcorrectCells);
        for j = 1 : ShuffledCrit.FcorrectNumberPassed
            cellNum = ShuffledCrit.FcorrectCells(j);
            figure(1)
            subplot(n,m,count)
            count = count +1;
            h = imagesc(passedcorrectF(:,:,j));
            axis('tight')
            colormap('parula')
            colorbar
            caxis([minmax(ShuffledCrit.FcorrectCells(j),1) minmax(ShuffledCrit.FcorrectCells(j),2)])
            view(2)
            xlabel('Time(frame)')
            ylabel('Trial Number')
            title(['Passed Criteria Correct: Cell ',num2str(cellNum)])
            if mod(j,n*m) == 0 || j == ShuffledCrit.FcorrectNumberPassed
                pause(0.01)
                count = 1;                
                saveas(gcf,[pwd,'/FrontAnchored/ChoiceSeperation/CorrectTrials/',num2str(choices(i)),'/', num2str(cellNum),'CorrectTrials.jpg']);
                clf
            end
        end
    else
        passedcorrectF = [];        
    end
    %Incorrect
    if ~isempty(SepFront.IncorrectTrial{choices(i)})
        passedincorrectF = SepFront.IncorrectTrial{choices(i)}(:,:,ShuffledCrit.FincorrectCells);
        for j = 1 : ShuffledCrit.FincorrectNumberPassed
            cellNum = ShuffledCrit.FincorrectCells(j);
            figure(1)
            subplot(n,m,count)
            count = count +1;
            h = imagesc(passedincorrectF(:,:,j));
            axis('tight')
            colormap('parula')
            colorbar
            caxis([minmax(ShuffledCrit.FincorrectCells(j),1) minmax(ShuffledCrit.FincorrectCells(j),2)])
            view(2)
            xlabel('Time(frame)')
            ylabel('Trial Number')
            title(['Passed Criteria Incorrect: Cell ',num2str(cellNum)])
            if mod(j,n*m) == 0 || j == ShuffledCrit.FincorrectNumberPassed
                pause(0.01)
                count = 1;
                saveas(gcf,[pwd,'/FrontAnchored/ChoiceSeperation/IncorrectTrials/',num2str(choices(i)),'/', num2str(cellNum),'IncorrectTrials.jpg']);
                clf
            end
        end
    else
        passedincorrectF = [];
    end
    %Correct Correction
    if ~isempty(SepFront.CorrectCorrectionTrial{choices(i)})
        passedcorrectcorrectionF = SepFront.CorrectCorrectionTrial{choices(i)}(:,:,ShuffledCrit.FccorCells);
        for j = 1 : ShuffledCrit.FccorNumberPassed
            cellNum = ShuffledCrit.FccorCells(j);
            figure(1)
            subplot(n,m,count)
            count = count +1;
            h = imagesc(passedcorrectcorrectionF(:,:,j));
            axis('tight')
            colormap('parula')
            colorbar
            caxis([minmax(ShuffledCrit.FccorCells(j),1) minmax(ShuffledCrit.FccorCells(j),2)])
            view(2)
            xlabel('Time(frame)')
            ylabel('Trial Number')
            title(['Passed Criteria Correct Correction: Cell ',num2str(cellNum)])
            if mod(j,n*m) == 0 || j == ShuffledCrit.FccorNumberPassed
                pause(0.01)
                count = 1;
                saveas(gcf,[pwd,'/FrontAnchored/ChoiceSeperation/CorrectCorrectionTrials/',num2str(choices(i)),'/', num2str(cellNum),'CorrectCorrectionTrials.jpg']);
                clf
            end
        end
    else
        passedcorrectcorrectionF = [];
    end
    %Incorrect Correction
    if IncCor && ~isempty(SepFront.IncorrectCorrectionTrial{choices(i)})
        passedincorrectcorrectionF = SepFront.IncorrectCorrectionTrial{choices(i)}(:,:,ShuffledCrit.FicorCells);
        for j = 1 : ShuffledCrit.FicorNumberPassed
            cellNum = ShuffledCrit.FicorCells(j);
            figure(1)
            subplot(n,m,count)
            count = count +1;
            h = imagesc(passedincorrectcorrectionF(:,:,j));
            axis('tight')
            colormap('parula')
            colorbar
            caxis([minmax(ShuffledCrit.FicorCells(j),1) minmax(ShuffledCrit.FicorCells(j),2)])
            view(2)
            xlabel('Time(frame)')
            ylabel('Trial Number')
            title(['Passed Criteria Incorrect Correction: Cell ',num2str(cellNum)])
            if mod(j,n*m) == 0 || j == ShuffledCrit.FicorNumberPassed
                pause(0.01)
                count = 1;
                saveas(gcf,[pwd,'/FrontAnchored/ChoiceSeperation/IncorrectCorrectionTrials/',num2str(choices(i)),'/', num2str(cellNum),'IncorrectCorrectionTrials.jpg']);
                clf
            end
        end
    else
        passedincorrectcorrectionF = [];
    end
    meancFpop = sortpeaks(permute(mean(passedcorrectF,1),[3 2 1]));
    meaniFpop = sortpeaks(permute(mean(passedincorrectF,1),[3 2 1]));
    meanccorFpop = sortpeaks(permute(mean(passedcorrectcorrectionF,1),[3 2 1]));
    if IncCor
        meanicorFpop = sortpeaks(permute(mean(passedincorrectcorrectionF,1),[3 2 1]));
    else
        meanicorFpop = [];
    end
    figure(1)
    subplot(2,2,1)
    imagesc(meancFpop)
    ylabel('Cell')
    xlabel('Time (frames)')
    title('Front Anchored: Correct Mean Fluorescence')
    subplot(2,2,2)
    imagesc(meaniFpop)
    ylabel('Cell')
    xlabel('Time (frames)')
    title('Front Anchored: Incorrect Mean Fluorescence')
    subplot(2,2,3)
    imagesc(meanccorFpop)
    ylabel('Cell')
    xlabel('Time (frames)')
    title('Front Anchored: Correct Correction Mean Fluorescence')
    subplot(2,2,4)
    imagesc(meanicorFpop)
    ylabel('Cell')
    xlabel('Time (frames)')
    title('Front Anchored: Incorrect Correction Mean Fluorescence')
    saveas(gcf,[pwd,'/FrontAnchored/ChoiceSeperation_',num2str(choices(i)),'_PopulationRasterPlot.jpg']);
    clf  
end

%Back Anchored
choices = PokeFind(SepBack);
count = 1;
for i = 1 : length(choices)
    %Correct
    if ~isempty(SepBack.CorrectTrial{choices(i)})
        passedcorrectB = SepBack.CorrectTrial{choices(i)}(:,:,ShuffledCrit.BcorrectCells);
        for j = 1 : ShuffledCrit.BcorrectNumberPassed
            cellNum = ShuffledCrit.BcorrectCells(j);
            figure(1)
            subplot(n,m,count)
            count = count +1;
            h = imagesc(passedcorrectB(:,:,j));
            axis('tight')
            colormap('parula')
            colorbar
            caxis([minmax(ShuffledCrit.BcorrectCells(j),1) minmax(ShuffledCrit.BcorrectCells(j),2)])
            view(2)
            xlabel('Time(frame)')
            ylabel('Trial Number')
            title(['Passed Criteria Correct: Cell ',num2str(cellNum)])
            if mod(j,n*m) == 0 || j == ShuffledCrit.BcorrectNumberPassed
                pause(0.01)
                count = 1;                
                saveas(gcf,[pwd,'/BackAnchored/ChoiceSeperation/CorrectTrials/',num2str(choices(i)),'/', num2str(cellNum),'CorrectTrials.jpg']);
                clf
            end
        end
    else
        passedcorrectB = [];
    end
    %Incorrect
    if ~isempty(SepBack.IncorrectTrial{choices(i)})
        passedincorrectB = SepBack.IncorrectTrial{choices(i)}(:,:,ShuffledCrit.BincorrectCells);
        for j = 1 : ShuffledCrit.BincorrectNumberPassed
            cellNum = ShuffledCrit.BincorrectCells(j);
            figure(1)
            subplot(n,m,count)
            count = count +1;
            h = imagesc(passedincorrectB(:,:,j));
            axis('tight')
            colormap('parula')
            colorbar
            caxis([minmax(ShuffledCrit.BincorrectCells(j),1) minmax(ShuffledCrit.BincorrectCells(j),2)])
            view(2)
            xlabel('Time(frame)')
            ylabel('Trial Number')
            title(['Passed Criteria Incorrect: Cell ',num2str(cellNum)])
            if mod(j,n*m) == 0 || j == ShuffledCrit.BincorrectNumberPassed
                pause(0.01)
                count = 1;
                saveas(gcf,[pwd,'/BackAnchored/ChoiceSeperation/IncorrectTrials/',num2str(choices(i)),'/', num2str(cellNum),'IncorrectTrials.jpg']);
                clf
            end
        end
    else
        passedincorrectB = [];
    end
    %Correct Correction
    if ~isempty(SepBack.CorrectCorrectionTrial{choices(i)})
        passedcorrectcorrectionB = SepBack.CorrectCorrectionTrial{choices(i)}(:,:,ShuffledCrit.BccorCells);
        for j = 1 : ShuffledCrit.BccorNumberPassed
            cellNum = ShuffledCrit.BccorCells(j);
            figure(1)
            subplot(n,m,count)
            count = count +1;
            h = imagesc(passedcorrectcorrectionB(:,:,j));
            axis('tight')
            colormap('parula')
            colorbar
            caxis([minmax(ShuffledCrit.BccorCells(j),1) minmax(ShuffledCrit.BccorCells(j),2)])
            view(2)
            xlabel('Time(frame)')
            ylabel('Trial Number')
            title(['Passed Criteria Correct Correction: Cell ',num2str(cellNum)])
            if mod(j,n*m) == 0 || j == ShuffledCrit.BccorNumberPassed
                pause(0.01)
                count = 1;
                saveas(gcf,[pwd,'/BackAnchored/ChoiceSeperation/CorrectCorrectionTrials/',num2str(choices(i)),'/', num2str(cellNum),'CorrectCorrectionTrials.jpg']);
                clf
            end
        end
    else
        passedcorrectcorrectionB = [];
    end
    %Incorrect Correction
    if IncCor && ~isempty(SepBack.IncorrectCorrectionTrial{choices(i)})
        passedincorrectcorrectionB = SepBack.IncorrectCorrectionTrial{choices(i)}(:,:,ShuffledCrit.BicorCells);
        for j = 1 : ShuffledCrit.BicorNumberPassed
            cellNum = ShuffledCrit.BicorCells(j);
            figure(1)
            subplot(n,m,count)
            count = count +1;
            h = imagesc(passedincorrectcorrectionB(:,:,j));
            axis('tight')
            colormap('parula') 
            colorbar
            caxis([minmax(ShuffledCrit.BicorCells(j),1) minmax(ShuffledCrit.BicorCells(j),2)])
            view(2)
            xlabel('Time(frame)')
            ylabel('Trial Number')
            title(['Passed Criteria Incorrect Correction: Cell ',num2str(cellNum)])
            if mod(j,n*m) == 0 || j == ShuffledCrit.BicorNumberPassed
                pause(0.01)
                count = 1;
                saveas(gcf,[pwd,'/BackAnchored/ChoiceSeperation/IncorrectCorrectionTrials/',num2str(choices(i)),'/', num2str(cellNum),'IncorrectCorrectionTrials.jpg']);
                clf
            end
        end
    else
        passedincorrectcorrectionB = [];
    end
    %Figure2: Sample selective activity in the population
    meancBpop = sortpeaks(permute(mean(passedcorrectB,1),[3 2 1]));
    meaniBpop = sortpeaks(permute(mean(passedincorrectB,1),[3 2 1]));
    meanccorBpop = sortpeaks(permute(mean(passedcorrectcorrectionB,1),[3 2 1]));
    if IncCor
        meanicorBpop = sortpeaks(permute(mean(passedincorrectcorrectionB,1),[3 2 1]));
    else
        meanicorBpop = [];
    end
    
    figure(1)
    subplot(2,2,1)
    imagesc(meancBpop)
    ylabel('Cell')
    xlabel('Time (frames)')
    title('Back Anchored: Correct Mean Fluorescence')
    subplot(2,2,2)
    imagesc(meaniBpop)
    ylabel('Cell')
    xlabel('Time (frames)')
    title('Back Anchored: Incorrect Mean Fluorescence')
    subplot(2,2,3)
    imagesc(meanccorBpop)
    ylabel('Cell')
    xlabel('Time (frames)')
    title('Back Anchored: Correct Correction Mean Fluorescence')
    subplot(2,2,4)
    imagesc(meanicorBpop)
    ylabel('Cell')
    xlabel('Time (frames)')
    title('Back Anchored: Incorrect Correction Mean Fluorescence')
    saveas(gcf,[pwd,'/BackAnchored/ChoiceSeperation_',num2str(choices(i)),'_PopulationRasterPlot.jpg']);
    clf  
end

%Step2: Correct vs Incorrect Trial comparison
%Figure1: Visual Comparison of correct vs all other Trials at the
%individual cell level
%Figure2: Visal comparison of correct vs all other tirals at the population
%level

%Delay
DisplayComp(SepDelay,minmax,'Delay',ShuffledCrit.DcorrectCells);
%Front Anchored
DisplayComp(SepFront,minmax,'FrontAnchored',ShuffledCrit.FcorrectCells);
%Back Anchored
DisplayComp(SepBack,minmax,'BackAnchored',ShuffledCrit.BcorrectCells);

clf
%Figure4: Locating the source of the error signal
%Delay
confInterD = ErrorSource(SepDelay,1,pwd,ShuffledCrit.DcorrectCells);
%Front Anchored
confInterF = ErrorSource(SepFront,2,pwd,ShuffledCrit.FcorrectCells);
%Back Anchored
confInterB = ErrorSource(SepBack,3,pwd,ShuffledCrit.BcorrectCells);

%Figure3: Comparison of correlation values between correct and incorrect
%trials

%Delay
metric1D = tMetric1(SepDelay,1);
metric2D = tMetric2(SepDelay);
count = 1;
for i = 1 : length(metric1D)
    for j = 1 : length(metric1D)
        if i == j && ~isempty(metric1D{i,j})
            figure(1)
            cdfplot(metric2D(:,count))
            hold on
            if ~isempty(metric1D{i,j}.Incorrect)
                cdfplot(metric1D{i,j}.Incorrect)
            end
            title('Delay CDF plot')
            legend('Metric2: Control Correct subset comparison','Metric1: correlation between incorrect and correct')
            count = count+1;
            saveas(gcf, [pwd,'/Delay/ConfidenceIntervals/CDFplot_Choice_',num2str(i),'.jpg'])
            pause(0.01);
            clf
        end
    end
end
%FrontAnchored
metric1F = tMetric1(SepFront,2);
metric2F = tMetric2(SepFront);
count = 1;
for i = 1 : length(metric1F)
    for j = 1 : length(metric1F)
        if i == j && ~isempty(metric1F{i,j})
            figure(1)
            cdfplot(metric2F(:,count))
            hold on
            if ~isempty(metric1F{i,j}.Incorrect)
                cdfplot(metric1F{i,j}.Incorrect)
            end
            
            title('Front Anchored CDF plot')
            legend('Metric2: Control Correct subset comparison','Metric1: correlation between incorrect and correct')
            count = count+1;
            saveas(gcf, [pwd,'/FrontAnchored/ConfidenceIntervals/CDFplot_Choice_',num2str(i),'.jpg'])
            pause(0.01);
            clf
        end
    end
end
%BackAnchored
metric1B = tMetric1(SepBack,3);
metric2B = tMetric2(SepBack);
count = 1;
for i = 1 : length(metric1B)
    for j = 1 : length(metric1B)
        if i == j && ~isempty(metric1B{i,j})
            figure(1)
            cdfplot(metric2B(:,count))
            hold on
            if ~isempty(metric1B{i,j}.Incorrect)
                cdfplot(metric1B{i,j}.Incorrect)
            end
            title('Back Anchored CDF plot')
            legend('Metric2: Control Correct subset comparison','Metric1: correlation between incorrect and correct')
            count = count+1;
            saveas(gcf, [pwd,'/BackAnchored/ConfidenceIntervals/CDFplot_Choice_',num2str(i),'.jpg'])
            pause(0.01);
            clf
        end
    end
end


%-------------------STEP 5------------------------

out = step_5(ms,behav,events,groupID,item_name,eventTime,eventInd); 

