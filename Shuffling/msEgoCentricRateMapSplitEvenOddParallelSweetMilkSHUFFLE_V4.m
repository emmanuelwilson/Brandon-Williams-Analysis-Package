function [out]= msEgoCentricRateMapSplitEvenOddParallelSweetMilkSHUFFLE_V4(ms,HD,tracking, frameMap, dimX, dimY, deconvolve, QPO, QPW,DistBinSize,DegBinSize)
%%Egocentric Boundary Cell Rate Map function,boundary location polar plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function will take your data and analize it in a way to facilitate  %
%Egocentric Boundary Cell(EBC) identification.                            %
%   INPUTS:                                                               %
%       - ms: Contains Calcium traces ms.deconvolvedSig and/ ms.FiltTraces%
%       - HD: Head tracking data, N x 1 matrix of head orientation        %
%       - Tracking: N x 2 matrix for X and Y position of mouse            %
%       - frameMap: Synchronisation matrix, matches behaviour frames with %
%       physiology frames.                                                %
%       -dimX: Dimension of open field in the X direction                 %
%       -dimY: Dimensions of open field in the Y direction                %
%       -deconvolve: 1 or 0 logical to use deconvoled signal over pure    %
%       calcium signal                                                    %
%       -QPO: coordinates of object corners, [] to manual set             %
%       -QPW: coordinates of Wall corners, [] to manual set               %
%   OUTPUT:                                                               %
%       - out: Contains all EgoCentric ratemaps, MRL values, analysis     %
%       times, correlation values and object/wall definitions.            %
%This function will create a new folder within you directory and save all %
%analysis figures and matrices in in said folder.                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author Emmanuel Wilson

name = 'EBCevenOddSplitParallelDeconvolved_Donut3SWEETMILK_D1d2A6_TEST_fix.mat';

%Set Analysis Distance
if dimX > dimY
    FOVsize = round(dimY/2);
else
    FOVsize = round(dimX/2);
end

%set deconvolved signal if being used
if deconvolve
    ms.FiltTraces = ms.deconvolvedSig;
end

%% Get behavior information
ratemapsO = zeros(FOVsize/DistBinSize,360/DegBinSize,length(ms.FiltTraces(1,:)));                  %Probability ratemap values
ratemaps1O = zeros(FOVsize/DistBinSize,360/DegBinSize,length(ms.FiltTraces(1,:)));                 %Probability ratemap values
ratemaps2O = zeros(FOVsize/DistBinSize,360/DegBinSize,length(ms.FiltTraces(1,:)));                 %Probability ratemap values

ratemapsW = zeros(FOVsize/DistBinSize,360/DegBinSize,length(ms.FiltTraces(1,:)));        %Probability ratemap values
ratemaps1W = zeros(FOVsize/DistBinSize,360/DegBinSize,length(ms.FiltTraces(1,:)));        %Probability ratemap values
ratemaps2W = zeros(FOVsize/DistBinSize,360/DegBinSize,length(ms.FiltTraces(1,:)));        %Probability ratemap values
mkdir(name)                                                                 %Create new folder within current directory
degBins = (-180:DegBinSize:179);                                                     %Angle bins for EBC metric polar plot
degBins = degBins';                                                         %reorient bins
degBins = deg2rad(degBins);                                                 %Convert to radians
distanceBins = 0:DistBinSize:FOVsize;                                                 %set to look at half the length of the field which in our case is ~38cm (37.5 rounded up)
fps = 30;                                                                   %Frames per second
spf = 1/fps;                                                                %Seconds per frame
ms.timestamp = frameMap.*spf;                                               %time stamp in seconds
minDist = [1:2];                                                            %minimum distance ignored

%% FrameMap special case correction
%IF SPECIAL CASE APPEARS MAKE SURE TO CHECK framemap, ms AND SINKdata ARE
%CORRECT
%If framemap exceeds number of frames available during tracking (incase of
%behavioural anomily where behav videos were recorded past the experiment)
if length(frameMap)> length(tracking(:,1))
    i = 0;
    test = frameMap(1:find(frameMap == length(tracking(:,1))));
    while isempty(test)
        test = frameMap(1:find(frameMap == length(tracking(:,1))-i));
        i = i +1;
    end
    frameMap = test;
    fprintf('SPECIAL CASE: FrameMap is larger than the behav')
    beep
    pause
end

%If miniscope recording used is longer than the trace being used (incase
%miniscope falls off or records past the useful experimental length)
if length(frameMap)> length(ms.FiltTraces(:,1))
    i = 0;
    test = frameMap(1:length(ms.FiltTraces(:,1)));
    while isempty(test)
        test = frameMap(1:length(ms.FiltTraces(:,1))-i);
        i = i +1;
    end
    frameMap = test;
    fprintf('SPECIAL CASE: FrameMap is larger than Physiology')
    beep
    pause
end

%% Get structure of environment
%Identify where the bounds of the environment are located. Through a subfunction that allows the user to click
%on a plot of the positional data to indicate where the corners of the environment are located.

%% Calculate distances
if isempty(QPO)
    QPO = findEdges(tracking);       %finds the edge of the graph and mark the boundaries
end

if isempty(QPW)
    QPW = findEdges(tracking);       %finds the edge of the graph and mark the boundaries
end

%Set pixel/cm ratio
pixX = dimX/(max(QPW(:,1)) - min(QPW(:,1)));
pixY = dimY/(max(QPW(:,2)) - min(QPW(:,2)));

%Set Figure size
set(groot, 'defaultFigureUnits','normalized');
set(groot, 'defaultFigurePosition',[0 0 1 1]);

%% Calculate distances
degSamp = 1;                                                                %angle resolution
[disW, ex, ey] = subfuncW(tracking(frameMap(:,1),1),tracking(frameMap(:,1),2),HD(frameMap), QPW, degSamp);   %calls funtion to bring back wall distances when neuron fired
dis_rawW = disW;
disW = fillmissing(disW,'pchip',2);                                         %interpolates missing values
disW = disW*pixX;                                                           %Converts boundary distances from pixels to cm.
disW = circshift(disW,90,2);                                                %shifts values by 90 degrees

degSamp = 1;                                                                %angle resolution
[disO, ex, ey] = subfuncO(tracking(frameMap(:,1),1),tracking(frameMap(:,1),2),HD(frameMap), QPO, degSamp);   %calls funtion to bring back wall distances when neuron fired
dis_rawO = disO;
disO = disO*pixX;                                                           %Converts boundary distances from pixels to cm.
disO = circshift(disO,90,2);               %shifts values by 90 degrees

firing = ms.deconvolvedSig;
parfor itteration = 1 : 100
    shuffledFiring = CShuffle(firing);                                  %Shuffle firing peaks
    distanceBinspar = distanceBins;
    mrlittO = zeros(length(firing(1,:)),1);
    mrlittO1 = zeros(length(firing(1,:)),1);
    mrlittO2 = zeros(length(firing(1,:)),1);
    mrlittW = zeros(length(firing(1,:)),1);
    mrlittW1 = zeros(length(firing(1,:)),1);
    mrlittW2 = zeros(length(firing(1,:)),1);
    %Loop through every cell, extract and analize firing instances and boundary locations
    for cellNum = 1 : length(ms.FiltTraces(1,:))                                                            %create parallel ms variable                                                         %Progress flag
        stime = tic;                                                            %processing time variable
        distanceBinsPar = distanceBins;                                         %parallel distance variable
        %Exctract calcium trace
        if deconvolve == 1
            fire = shuffledFiring(:,cellNum);
        else
            fire = shuffledFiring(:,cellNum);
        end
        fire = fire - min(fire);
        
        %Seperate Even and Odd minutes
        mins = length(fire)/1800;
        if mins <= round(mins)
            mins = round(mins);
        else
            mins = round(mins)+1;
        end
        for i = 1 : mins
            if mod(i,2) == 0
                if i == 2
                    fire2 = fire((1800*(i-1)+1:1800*i),1);
                    iTime2 = (1800*(i-1)+1:1800*i);
                    dis2 = disW((1800*(i-1)+1:1800*i),:);
                    dis2O = disO((1800*(i-1)+1:1800*i),:);
                elseif i < mins
                    fire2(end+1:end + 1800) = fire(1800*(i-1)+1:1800*i,1);
                    iTime2(end+1:end + 1800) =(1800*(i-1)+1:1800*i);
                    dis2(end+1:end + 1800,:) = disW(1800*(i-1)+1:1800*i,:);
                    dis2O(end+1:end + 1800,:) = disO(1800*(i-1)+1:1800*i,:);
                else
                    fire2(end+1:end+(length(fire)-((mins-1)*1800))) = fire((1800*(i-1)+1:end),1);
                    iTime2(end+1:end+(length(fire)-((mins-1)*1800))) = (1800*(i-1)+1:length(fire));
                    dis2(end+1:end+(length(fire)-((mins-1)*1800)),:) = disW((1800*(i-1)+1:end),:);
                    dis2O(end+1:end+(length(fire)-((mins-1)*1800)),:) = disO((1800*(i-1)+1:end),:);
                end
            else
                if i == 1
                    fire1 = fire(1:1800*i,1);
                    iTime1 = (1:1800*i);
                    dis1 = disW(1:1800,:);
                    dis1O = disO(1:1800,:);
                elseif i < mins
                    fire1(end+1:end+1800) = fire(((1800*i)-1799:1800*i),1);
                    iTime1(end+1:end+1800)= ((1800*i)-1799:1800*i);
                    dis1(end+1:end+1800,:) = disW((1800*i)-1799:1800*i,:);
                    dis1O(end+1:end+1800,:) = disO((1800*i)-1799:1800*i,:);
                else
                    fire1(end+1:end+(length(fire)-(mins-1)*1800)) = fire(((1800*i)-1799:end),1);
                    iTime1(end+1:end+(length(fire)-(mins-1)*1800)) = ((1800*i)-1799:length(fire));
                    dis1(end+1:end+(length(fire)-(mins-1)*1800),:) = disW((1800*i)-1799:length(fire),:);
                    dis1O(end+1:end+(length(fire)-(mins-1)*1800),:) = disO((1800*i)-1799:length(fire),:);
                end
            end
        end
        
        ifire = find(fire);                                                     %Find indices for all non-zero values
        ifire1 = find(fire1);                                                   %Find indices for all non-zero values
        ifire2 = find(fire2);                                                   %Find indices for all non-zero values
        
        if(~isempty(ifire) && ~isempty(ifire1) && ~isempty(ifire2))
            for j = 1 : length(ifire)
                %full run
                ind_fire(j) = ifire(j);                                         %firing index
                %first half
                if j < length(ifire1)
                    ind_fire1(j) = ifire1(j);                                   %firing index
                end
                %second half
                if j < length(ifire2)
                    ind_fire2(j) = ifire2(j);                                   %firing index
                end
            end
            
            %% Calculate raw maps:
            thetaBins = deg2rad(linspace(-180,180,size(disW,2)));               %angle bins
            thetaBins3d = deg2rad(linspace(-180,180,round(size(disW,2)/DegBinSize)));                    %angle bins
            occO = NaN(length(thetaBins), length(distanceBinsPar));                %Object total occupancy bins
            nspkO = occO;                                                       %Object occupancy during activity bins
            nspk1O = nspkO;
            nspk2O = nspkO;
            occ1O = occO;
            occ2O = occO;
            
            occW = NaN(length(thetaBins), length(distanceBinsPar));                %wall occupancy bins
            nspkW = occW;                                                       %Wall occupancy during activity bins
            nspk1W = nspkW;
            nspk2W = nspkW;
            occ1W = occW;
            occ2W = occW;
            
            distanceBinsPar(end+1) = Inf;                                       %Adds an Infinity value at the end of the bins as safety procaution/break point
            ci = ind_fire(:);                                                   %firing instances of the cell total length
            ci1 = ind_fire1(:);                                                 %firing instances of the cell odd mins
            ci2 = ind_fire2(:);                                                 %firing instances of the cell even mins
            for i = 1:length(thetaBins)
                t = disW(:,i);                                                  %boundary distance for a particular bin
                t1 = dis1(:,i);
                t2 = dis2(:,i);
                
                tO = disO(:,i); %boundary distance for a particular bin
                tO1 = dis1O(:,i);
                tO2 = dis2O(:,i);
                
                for k = 1:length(distanceBinsPar)-1
                    %Object
                    inds = tO>=distanceBinsPar(k) & tO<distanceBinsPar(k+1);    %filter through the boundary distances
                    occO(i,k) = sum(inds);                                      %Object occupancy definition
                    inds = find(inds);                                          %find all non-zero boundary distances indices
                    nspkO(i,k) = sum(fire(intersect(inds,ci)));                 %Number of spike instances definition
                    %first half
                    inds1 = tO1>=distanceBinsPar(k) & tO1<distanceBinsPar(k+1); %filter through the boundary distances
                    occ1O(i,k) = sum(inds1);                                    %Object occupancy definition
                    inds1 = find(inds1);                                        %find all non-zero boundary distances indices
                    nspk1O(i,k) = sum(fire1(intersect(inds1,ci1)));             %Number of spike instances definition
                    %second half
                    inds2 = tO2>=distanceBinsPar(k) & tO2<distanceBinsPar(k+1); %filter through the boundary distances
                    occ2O(i,k) = sum(inds2);                                    %Object occupancy definition
                    inds2 = find(inds2);                                        %find all non-zero boundary distances indices
                    inds2 = inds2;
                    nspk2O(i,k) = sum(fire2(intersect(inds2,ci2)));             %Number of spike instances definition
                    
                    %Wall
                    inds = t>=distanceBinsPar(k) & t<distanceBinsPar(k+1);      %filter through the boundary distances
                    occW(i,k) = sum(inds);                                      %Wall occupancy definition
                    inds = find(inds);                                          %find all non-zero boundary distances indices
                    nspkW(i,k) = sum(fire(intersect(inds,ci)));                 %Number of spike instances definition
                    %first half
                    inds1 = t1>=distanceBinsPar(k) & t1<distanceBinsPar(k+1);   %filter through the boundary distances
                    occ1W(i,k) = sum(inds1);                                    %Wall occupancy definition
                    inds1 = find(inds1);                                        %find all non-zero boundary distances indices
                    nspk1W(i,k) = sum(fire1(intersect(inds1,ci1)));             %Number of spike instances definition
                    %second half
                    inds2 = t2>=distanceBinsPar(k) & t2<distanceBinsPar(k+1);   %filter through the boundary distances
                    occ2W(i,k) = sum(inds2);                                    %Wall occupancy definition
                    inds2 = find(inds2);                                        %find all non-zero boundary distances indices
                    inds2 = inds2;
                    nspk2W(i,k) = sum(fire2(intersect(inds2,ci2)));             %Number of spike instances definition
                end
            end
            %convert from 1 degree to 3 degree bins
            occ3dO = zeros(360/DegBinSize,length(occO(1,:)));
            occ13dO = zeros(360/DegBinSize,length(occO(1,:)));
            occ23dO = zeros(360/DegBinSize,length(occO(1,:)));
            
            nspk3dO = zeros(360/DegBinSize,length(occO(1,:)));
            nspk13dO = zeros(360/DegBinSize,length(occO(1,:)));
            nspk23dO = zeros(360/DegBinSize,length(occO(1,:)));
            
            occ3dW = zeros(360/DegBinSize,length(occW(1,:)));
            occ13dW = zeros(360/DegBinSize,length(occW(1,:)));
            occ23dW = zeros(360/DegBinSize,length(occW(1,:)));
            
            nspk3dW = zeros(360/DegBinSize,length(occW(1,:)));
            nspk13dW = zeros(360/DegBinSize,length(occW(1,:)));
            nspk23dW = zeros(360/DegBinSize,length(occW(1,:)));
            
            for i = 1 : round(length(thetaBins)/DegBinSize)
                octempO = sum(occO(i*DegBinSize-(DegBinSize-1):i*DegBinSize,:),1);
                octemp1O = sum(occ1O(i*DegBinSize-(DegBinSize-1):i*DegBinSize,:),1);
                octemp2O = sum(occ2O(i*DegBinSize-(DegBinSize-1):i*DegBinSize,:),1);
                
                occ3dO(i,:) = octempO;
                occ13dO(i,:) = octemp1O;
                occ23dO(i,:) = octemp2O;
                
                nspktempO = sum(nspkO(i*DegBinSize-(DegBinSize-1):i*DegBinSize,:),1);
                nspktemp1O= sum(nspk1O(i*DegBinSize-(DegBinSize-1):i*DegBinSize,:),1);
                nspktemp2O = sum(nspk2O(i*DegBinSize-(DegBinSize-1):i*DegBinSize,:),1);
                
                nspk3dO(i,:) = nspktempO;
                nspk13dO(i,:) = nspktemp1O;
                nspk23dO(i,:) = nspktemp2O;
                
                octempW = sum(occW(i*DegBinSize-(DegBinSize-1):i*DegBinSize,:),1);
                octemp1W = sum(occ1W(i*DegBinSize-(DegBinSize-1):i*DegBinSize,:),1);
                octemp2W = sum(occ2W(i*DegBinSize-(DegBinSize-1):i*DegBinSize,:),1);
                
                occ3dW(i,:) = octempW;
                occ13dW(i,:) = octemp1W;
                occ23dW(i,:) = octemp2W;
                
                nspktempW = sum(nspkW(i*DegBinSize-(DegBinSize-1):i*DegBinSize,:),1);
                nspktemp1W= sum(nspk1W(i*DegBinSize-(DegBinSize-1):i*DegBinSize,:),1);
                nspktemp2W = sum(nspk2W(i*DegBinSize-(DegBinSize-1):i*DegBinSize,:),1);
                
                nspk3dW(i,:) = nspktempW;
                nspk13dW(i,:) = nspktemp1W;
                nspk23dW(i,:) = nspktemp2W;
            end
            
            occO = occ3dO;
            occ1O = occ13dO;
            occ2O = occ23dO;
            
            nspkO = nspk3dO;
            nspk1O = nspk13dO;
            nspk2O = nspk23dO;
            
            occW = occ3dW;
            occ1W = occ13dW;
            occ2W = occ23dW;
            
            nspkW = nspk3dW;
            nspk1W = nspk13dW;
            nspk2W = nspk23dW;
            
            distanceBinsPar = distanceBinsPar(1:end-1);                         %bring back distance bins to not go to infinity
            if any(nspkO(:)>0) && any(nspkW(:)>0)
                %Object
                % bring back to original dims
                occO=occO(:,1:end-1);occO=occO';
                cutout = find(occO<0);
                occO(cutout) = 0;
                nspkO=nspkO(:,1:end-1);nspkO=nspkO';
                %first half
                occ1O=occ1O(:,1:end-1);occ1O=occ1O';
                occ1O(cutout) = 0;
                nspk1O=nspk1O(:,1:end-1);nspk1O=nspk1O';
                %second half
                occ2O=occ2O(:,1:end-1);occ2O=occ2O';
                occ2O(cutout) = 0;
                nspk2O=nspk2O(:,1:end-1);nspk2O=nspk2O';
                
                rmO = (nspkO./occO) * fps;
                rmO(find(isnan(rmO))) = min(rmO(:));
                rmO(find(isinf(rmO))) = min(rmO(:));
                rmO = rmO - min(rmO(:));
                
                rm1O = (nspk1O./occ1O) * fps;
                rm1O(find(isnan(rm1O))) = min(rm1O(:));
                rm1O(find(isinf(rm1O))) = min(rm1O(:));
                rm1O = rm1O - min(rm1O(:));
                
                rm2O = (nspk2O./occ2O) * fps;
                rm2O(find(isnan(rm2O))) = min(rm2O(:));
                rm2O(find(isinf(rm2O))) = min(rm2O(:));
                rm2O = rm2O - min(rm2O(:));
                
                %Wall
                occW=occW(:,1:end-1);occW=occW';
                cutout = find(occW<0);
                occW(cutout) = 0;
                nspkW=nspkW(:,1:end-1);nspkW=nspkW';
                %first half
                occ1W=occ1W(:,1:end-1);occ1W=occ1W';
                occ1W(cutout) = 0;
                nspk1W=nspk1W(:,1:end-1);nspk1W=nspk1W';
                %second half
                occ2W=occ2W(:,1:end-1);occ2W=occ2W';
                occ2W(cutout) = 0;
                nspk2W=nspk2W(:,1:end-1);nspk2W=nspk2W';
                
                rmW = (nspkW./occW) * fps;
                rmW(find(isnan(rmW))) = min(rmW(:));
                rmW(find(isinf(rmW))) = min(rmW(:));
                rmW = rmW - min(rmW(:));
                
                rm1W = (nspk1W./occ1W) * fps;
                rm1W(find(isnan(rm1W))) = min(rm1W(:));
                rm1W(find(isinf(rm1W))) = min(rm1W(:));
                rm1W = rm1W - min(rm1W(:));
                
                rm2W = (nspk2W./occ2W) * fps;
                rm2W(find(isnan(rm2W))) = min(rm2W(:));
                rm2W(find(isinf(rm2W))) = min(rm2W(:));
                rm2W = rm2W - min(rm2W(:));
                
                %% Smoothing
                %ratemap
                %OBJECT
                %full run
                nd = numel(thetaBins3d);
                rmO = [rmO rmO rmO];
                rmO = CMBHOME.Utils.SmoothMat(rmO,smooth(1:2),smooth(3));   % Smooth it
                rmO = rmO(:,nd+1:2*nd); % bring it back
                %first half
                rm1O = [rm1O rm1O rm1O];
                rm1O = CMBHOME.Utils.SmoothMat(rm1O,smooth(1:2),smooth(3));   % Smooth it
                rm1O = rm1O(:,nd+1:2*nd); % bring it back
                %second half
                rm2O = [rm2O rm2O rm2O];
                rm2O = CMBHOME.Utils.SmoothMat(rm2O,smooth(1:2),smooth(3));   % Smooth it
                rm2O = rm2O(:,nd+1:2*nd); % bring it back
                
                rmO = fliplr(rmO);
                rmO(minDist,:) = 0;
                rm1O = fliplr(rm1O);
                rm1O(minDist,:) = 0;
                rm2O = fliplr(rm2O);
                rm2O(minDist,:) = 0;
                
                %WALL
                nd = numel(thetaBins3d);
                rmW = [rmW rmW rmW];
                rmW = CMBHOME.Utils.SmoothMat(rmW,smooth(1:2),smooth(3));   % Smooth it
                rmW = rmW(:,nd+1:2*nd); % bring it back
                %first half
                rm1W = [rm1W rm1W rm1W];
                rm1W = CMBHOME.Utils.SmoothMat(rm1W,smooth(1:2),smooth(3));   % Smooth it
                rm1W = rm1W(:,nd+1:2*nd); % bring it back
                %second half
                rm2W = [rm2W rm2W rm2W];
                rm2W = CMBHOME.Utils.SmoothMat(rm2W,smooth(1:2),smooth(3));   % Smooth it
                rm2W = rm2W(:,nd+1:2*nd); % bring it back
                
                rmW = fliplr(rmW);
                rmW(minDist,:) = 0;
                rm1W = fliplr(rm1W);
                rm1W(minDist,:) = 0;
                rm2W = fliplr(rm2W);
                rm2W(minDist,:) = 0;
                
                %% EBC METRIC
                %Object EBC Metric
                metricO = mean(rmO,1)';
                if ~deconvolve
                    metricO = metricO - min(metricO);
                end
                
                xs = metricO(1:end-1).*cos(degBins(1:end-1)); % average
                ys = metricO(1:end-1).*sin(degBins(1:end-1));
                
                ang_hdO = atan2(mean(ys),mean(xs)); % mean direction
                
                mrO = (cos(ang_hdO)*sum(xs) + sin(ang_hdO)*sum(ys)) / sum(metricO(1:end-1)); % mean resultant length
                
                %Wall EBC Metric
                metricW = mean(rmW,1)';
                if ~deconvolve
                    metricW = metricW - min(metricW);
                end
                
                xs = metricW(1:end-1).*cos(degBins(1:end-1)); % average
                ys = metricW(1:end-1).*sin(degBins(1:end-1));
                
                ang_hdW = atan2(mean(ys),mean(xs)); % mean direction
                
                mrW = (cos(ang_hdW)*sum(xs) + sin(ang_hdW)*sum(ys)) / sum(metricW(1:end-1)); % mean resultant length
                
                %Object EBC Metric
                metric1O = mean(rm1O,1)';
                if ~deconvolve
                    metric1O = metric1O - min(metric1O);
                end
                
                xs = metric1O(1:end-1).*cos(degBins(1:end-1)); % average
                ys = metric1O(1:end-1).*sin(degBins(1:end-1));
                
                ang_hd1O = atan2(mean(ys),mean(xs)); % mean direction
                
                mr1O = (cos(ang_hd1O)*sum(xs) + sin(ang_hd1O)*sum(ys)) / sum(metric1O(1:end-1)); % mean resultant length
                
                %Wall EBC Metric
                metric1W = mean(rm1W,1)';
                if ~deconvolve
                    metric1W = metric1W - min(metric1W);
                end
                
                xs = metric1W(1:end-1).*cos(degBins(1:end-1)); % average
                ys = metric1W(1:end-1).*sin(degBins(1:end-1));
                
                ang_hd1W = atan2(mean(ys),mean(xs)); % mean direction
                
                mr1W = (cos(ang_hd1W)*sum(xs) + sin(ang_hd1W)*sum(ys)) / sum(metric1W(1:end-1)); % mean resultant length
                
                %Object EBC Metric
                metric2O = mean(rm2O,1)';
                if ~deconvolve
                    metric2O = metric2O - min(metric2O);
                end
                
                xs = metric2O(1:end-1).*cos(degBins(1:end-1)); % average
                ys = metric2O(1:end-1).*sin(degBins(1:end-1));
                
                ang_hd2O = atan2(mean(ys),mean(xs)); % mean direction
                
                mr2O = (cos(ang_hd2O)*sum(xs) + sin(ang_hd2O)*sum(ys)) / sum(metric2O(1:end-1)); % mean resultant length
                
                %Wall EBC Metric
                metric2W = mean(rm2W,1)';
                if ~deconvolve
                    metric2W = metric2W - min(metric2W);
                end
                
                xs = metric2W(1:end-1).*cos(degBins(1:end-1)); % average
                ys = metric2W(1:end-1).*sin(degBins(1:end-1));
                
                ang_hd2W = atan2(mean(ys),mean(xs)); % mean direction
                
                mr2W = (cos(ang_hd2W)*sum(xs) + sin(ang_hd2W)*sum(ys)) / sum(metric2W(1:end-1)); % mean resultant length
                
                %Save Results
                
                mrlittO(cellNum) = mrO; %save the MRL
                mrlittO1(cellNum) =mr1O;
                mrlittO2(cellNum) =mr2O;
                
                mrlittW(cellNum) = mrW; %save the MRL
                mrlittW1(cellNum) =mr1W;
                mrlittW2(cellNum) =mr2W;                
                                
                ind_fire = NaN(length(shuffledFiring(:,1)),length(shuffledFiring(1,:)));        %Indices of neuron activity/firing                                                            %Indices of activity for head direction
                ind_fire1 = NaN(size(ind_fire1));                                
                ind_fire2 = NaN(size(ind_fire2));                
            end
        end
    end
    MRLsaveO(:,itteration) = mrlittO;
    MRLsaveO1(:,itteration)= mrlittO1;
    MRLsaveO2(:,itteration) = mrlittO2;
    MRLsaveW(:,itteration) = mrlittW;
    MRLsaveW1(:,itteration)= mrlittW1;
    MRLsaveW2(:,itteration) = mrlittW2;
end
%output Variable
out.mrallO = MRLsaveO(:);
out.mrallO1 = MRLsaveO1(:);
out.mrallO2 = MRLsaveO2(:);
out.percentil99thO = prctile(out.mrallO,99);
out.mrallW = MRLsaveW(:);
out.mrallW1 = MRLsaveW1(:);
out.mrallW2 = MRLsaveW2(:);
out.percentil99thW = prctile(out.mrallW,99);
out.dimX = dimX;
out.dimY = dimY;
out.QPO = QPO;
out.QPW = QPW;
out.dir = pwd;

save([pwd,'/','ObjectShuffledStats.mat'],'out');
end

%% Subfunctions

%This function calculates the distance from the animal to boundaries of the environment at each behavioral data point.
%The distance calculation has to be done for all orientations around the animal centered on the animal’s
%current heading direction. That is to say that the animal’s current heading is always 0° and the distance
%to the boundaries is calculated for each of the 360 one-degree bins around the animal.

function [dis, ex, ey] = subfuncO(rx,ry,hd, QP, degSamp)

mxd = sqrt((max(rx)-min(rx))^2 + (max(ry)-min(ry))^2);
degs = deg2rad(-180:degSamp:180);
hd = deg2rad(hd);

edg = splitter(QP);
edg = cell2mat(edg(:));
dis = NaN(numel(rx),size(edg,1), numel(degs));
dir = dis;

for i = 1:size(edg,1)
    x1=edg(i,1,1);x2=edg(i,1,2);
    y1=edg(i,2,1);y2=edg(i,2,2);
    
    for h = 1:numel(degs)
        hdof=degs(h);
        y3=ry;x3=rx;
        y4=ry+mxd*sin(hd+hdof);
        x4=rx+mxd*cos(hd+hdof);
        
        %https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection#Intersection_of_two_lines
        px1 = (x1.*y2-y1.*x2).*(x3-x4) - (x1-x2).*(x3.*y4-y3.*x4);
        px2 = (x1-x2).*(y3-y4) - (y1-y2).*(x3-x4);
        px  = px1./px2;
        
        py1 = (x1.*y2-y1.*x2).*(y3-y4) - (y1-y2).*(x3.*y4-y3.*x4);
        py2 = (x1-x2).*(y3-y4) - (y1-y2).*(x3-x4);
        py = py1./py2;
        
        d = sqrt((ry-py).^2 + (rx-px).^2);
        dis(:,i,h) = d;
        
        % need to filter down to the right direction ...
        dir(:,i,h) = wrapToPi(atan2(py-ry,px-rx)-(hd+hdof));
        
        % oh ... we were allowing forever.... filter by bounding box
        bb = [min(QP(:,1)) max(QP(:,1)); min(QP(:,2)) max(QP(:,2))];
        % |xmin, xmax|
        % |ymin, ymax|
        bb = round(bb,10);
        py = round(py,10);
        px = round(px,10);
        %         indexes = (px>bb(1,1) & px<bb(1,2) & py>bb(2,1) & py<bb(2,2));
        indexes = ~((px >= bb(1,1) & px <= bb(1,2)) & (py >= bb(2,1) & py <= bb(2,2)));
        dis(indexes,i,h) = NaN;
    end
    
end

dis(dis>mxd) = NaN;
dis(abs(dir)>pi/4) = NaN;

%% output
dis=squeeze(nanmin(dis,[],2));
for i = 1 :length(rx)
    if(rx(i)<max(edg(:,1,1)) && rx(i)>min(edg(:,1,1)) && ry(i)<max(edg(:,2,1)) && ry(i)>min(edg(:,2,1)))
        dis(i,:) = NaN;
    end
end
dd=repmat(degs,size(rx,1),1) + repmat(hd,1,numel(degs));
dx=dis.*cos(dd); dy=dis.*sin(dd);
ey=dy+repmat(ry,1,numel(degs));
ex=dx+repmat(rx,1,numel(degs));

end

function [dis, ex, ey] = subfuncW(rx,ry,hd, QP, degSamp)

mxd = sqrt((max(rx)-min(rx))^2 + (max(ry)-min(ry))^2);                  %sets bin radial maximum
degs = deg2rad(-180:degSamp:180);
hd = deg2rad(hd);

edg = splitter(QP);
edg = cell2mat(edg(:));
dis = NaN(numel(rx),size(edg,1), numel(degs));
dir = dis;

for i = 1:size(edg,1)
    x1=edg(i,1,1);x2=edg(i,1,2);
    y1=edg(i,2,1);y2=edg(i,2,2);
    
    for h = 1:numel(degs)
        hdof=degs(h);
        y3=ry;x3=rx;
        y4=ry+mxd*sin(hd+hdof);
        x4=rx+mxd*cos(hd+hdof);
        
        %https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection#Intersection_of_two_lines
        px1 = (x1.*y2-y1.*x2).*(x3-x4) - (x1-x2).*(x3.*y4-y3.*x4);
        px2 = (x1-x2).*(y3-y4) - (y1-y2).*(x3-x4);
        px  = px1./px2;
        
        py1 = (x1.*y2-y1.*x2).*(y3-y4) - (y1-y2).*(x3.*y4-y3.*x4);
        py2 = (x1-x2).*(y3-y4) - (y1-y2).*(x3-x4);
        py = py1./py2;
        
        d = sqrt((ry-py).^2 + (rx-px).^2);
        dis(:,i,h) = d;
        
        % need to filter down to the right direction ...
        dir(:,i,h) = wrapToPi(atan2(py-ry,px-rx)-(hd+hdof));
        
        % oh ... we were allowing forever.... filter by bounding box
        bb = [min(QP(:,1)) max(QP(:,1)); min(QP(:,2)) max(QP(:,2))];
        % |xmin, xmax|
        % |ymin, ymax|
        indexes = ~(px>=bb(1,1) & px<=bb(1,2) & py>=bb(2,1) & py<=bb(2,2));
        dis(indexes,i,h) = NaN;
    end
    
end


dis(dis>mxd) = NaN;
dis(abs(dir)>pi/4) = NaN;

%% output
dis=squeeze(nanmin(dis,[],2));
for i = 1 :length(rx)
    if(rx(i)>max(edg(:,1,1)) || rx(i)<min(edg(:,1,1)) || ry(i)>max(edg(:,2,1)) || ry(i)<min(edg(:,2,1)))
        dis(i,:) = NaN;
    end
end
dd=repmat(degs,size(rx,1),1) + repmat(hd,1,numel(degs));
dx=dis.*cos(dd); dy=dis.*sin(dd);
ey=dy+repmat(ry,1,numel(degs));
ex=dx+repmat(rx,1,numel(degs));

end

%This subfunction will ask for the corner locations to determine the open
%field
function QP = findEdges(tracking)
ifEscape = 0;
h=figure();

while ~ifEscape
    figure(h);
    clf
    
    %[occupancy, xdim, ydim]=root.Occupancy([],[],1,2);
    %imagesc(xdim,ydim,occupancy');
    set(gca,'YDir','Normal'); %colormap(jet);
    clim=get(gca,'clim');set(gca,'clim',clim/50);
    hold on
    plot(tracking(:,1),tracking(:,2),'k');
    QP = [];
    
    set(h,'Name','Select Corners of Walls. Esc--> done. **Do not complete!**')
    
    button = 1;
    
    while button~=27
        [x,y,button] = ginput(1);
        
        clf
        
        %imagesc(xdim,ydim,occupancy');
        set(gca,'YDir','Normal'); %colormap(jet);
        clim=get(gca,'clim');set(gca,'clim',clim/50);
        hold on
        plot(tracking(:,1),tracking(:,2),'k');
        
        if ~isempty(QP)
            plot(QP(:,1),QP(:,2),'r')
            plot(QP(:,1),QP(:,2),'ro','MarkerFaceColor','r')
        end
        
        if button == 32 %space bar
            QP = [QP; NaN NaN];
        elseif button~=27
            QP = [QP; x y];
        end
        
        plot(QP(:,1),QP(:,2),'r')
        plot(QP(:,1),QP(:,2),'ro','MarkerFaceColor','r')
        
    end
    
    %Ask for verification
    edg = splitter(QP);
    clf;
    set(h,'Name','Verify. 0--> Try again; 1--> Confirm')
    plot(tracking(:,1),tracking(:,2),'k');
    hold on
    
    for m = 1:numel(edg)
        for n = 1:size(edg{m},1)
            sp = squeeze(edg{m}(n,:,1));
            ep = squeeze(edg{m}(n,:,2));
            plot([sp(1) ep(1)],[sp(2) ep(2)],'ro','MarkerFaceColor','r')
            plot([sp(1) ep(1)],[sp(2) ep(2)],'r')
        end
    end
    
    
    % set or repeat
    while button ~=48 && button~=49
        [~,~,button]=ginput(1);
    end
    ifEscape = button==49;
    
end

close(h);
drawnow();
end

%Split the corner coordinates in X and Y vectors
function edg = splitter(QP)

inds = find(isnan(QP(:,1)));
xs=CMBHOME.Utils.SplitVec(QP(:,1), @(x) isnan(x));
ys=CMBHOME.Utils.SplitVec(QP(:,2), @(x) isnan(x));

% split corners
for m = 1:size(xs,1)
    QP2{m} = [xs{m} ys{m}];
    QP2{m}(find(isnan(QP2{m}(:,1))),:) = [];
end

for m = 1:numel(QP2)
    for n = 1:size(QP2{m},1)
        sp = n;ep=n+1;
        if ep>size(QP2{m},1), ep=1;end
        edg{m}(n,:,1) = [QP2{m}(sp,1) QP2{m}(sp,2)];
        edg{m}(n,:,2) = [QP2{m}(ep,1) QP2{m}(ep,2)];
    end
end

end