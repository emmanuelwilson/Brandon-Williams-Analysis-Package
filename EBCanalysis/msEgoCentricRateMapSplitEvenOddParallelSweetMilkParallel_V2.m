function [out]= msEgoCentricRateMapSplitEvenOddParallelSweetMilkParallel_V2(ms,HD1,HD2,tracking1,tracking2, frameMap1,frameMap2, dimX, dimY, deconvolve, QPOL,QPOR, QPW,DistBinSize,DegBinSize,badframe1,badframe2)
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

name = 'Trial1Trial2ObjectData_D1A3_bothOb_V2';

%Set Analysis Distance 
if dimX > dimY
    FOVsize = round(dimY/2);
else
    FOVsize = round(dimX/2);
end

if ~isempty(badframe1) && ~isempty(badframe2)
    ms.FiltTraces(length(frameMap1):length(frameMap1)+ badframe2,:) = [];
    ms.FiltTraces(1:badframe1,:) = [];
    ms.deconvolvedSig(length(frameMap1):length(frameMap1)+ badframe2,:) = [];
    ms.deconvolvedSig(1:badframe1,:) = [];
    frameMap1(1:badframe1) = [];
    frameMap2(1:badframe2) = [];
elseif ~isempty(badframe1)
    ms.FiltTraces(1:badframe1,:) = [];
    ms.deconvolvedSig(1:badframe1,:) = [];
    frameMap1(1:badframe1) = [];
elseif ~isempty(badframe2)
    ms.FiltTraces(length(frameMap1):length(frameMap1)+ badframe2,:) = [];
    ms.deconvolvedSig(length(frameMap1):length(frameMap1)+ badframe2,:) = [];
    frameMap2(1:badframe2) = [];
end

%set deconvolved signal if being used
if deconvolve
    ms.FiltTraces = ms.deconvolvedSig;
end

%% Get behavior information

ratemapsOL = zeros(FOVsize/DistBinSize,360/DegBinSize,length(ms.FiltTraces(1,:)));                  %Probability ratemap values
ratemaps1OL = zeros(FOVsize/DistBinSize,360/DegBinSize,length(ms.FiltTraces(1,:)));                 %Probability ratemap values
ratemaps2OL = zeros(FOVsize/DistBinSize,360/DegBinSize,length(ms.FiltTraces(1,:)));                 %Probability ratemap values

ratemapsOR = zeros(FOVsize/DistBinSize,360/DegBinSize,length(ms.FiltTraces(1,:)));                  %Probability ratemap values
ratemaps1OR = zeros(FOVsize/DistBinSize,360/DegBinSize,length(ms.FiltTraces(1,:)));                 %Probability ratemap values
ratemaps2OR = zeros(FOVsize/DistBinSize,360/DegBinSize,length(ms.FiltTraces(1,:)));                 %Probability ratemap values

ratemapsW = zeros(FOVsize/DistBinSize,360/DegBinSize,length(ms.FiltTraces(1,:)));        %Probability ratemap values
ratemaps1W = zeros(FOVsize/DistBinSize,360/DegBinSize,length(ms.FiltTraces(1,:)));        %Probability ratemap values
ratemaps2W = zeros(FOVsize/DistBinSize,360/DegBinSize,length(ms.FiltTraces(1,:)));        %Probability ratemap values
mkdir(name)                                                                 %Create new folder within current directory
mkdir([name '\wall'])
mkdir([name '\ObjectSeperate'])
mkdir([name '\ObjectsTogether'])
degBins = (-180:DegBinSize:179);                                                     %Angle bins for EBC metric polar plot
degBins = degBins';                                                         %reorient bins
degBins = deg2rad(degBins);                                                 %Convert to radians
distanceBins = 0:DistBinSize:FOVsize;                                                 %set to look at half the length of the field which in our case is ~38cm (37.5 rounded up)
fps = 30;                                                                   %Frames per second
spf = 1/fps;                                                                %Seconds per frame
% ms.timestamp = fram eMap.*spf;                                               %time stamp in seconds
minDist = [1];                                                            %minimum distance ignored

%% FrameMap special case correction
%IF SPECIAL CASE APPEARS MAKE SURE TO CHECK framemap, ms AND SINKdata ARE
%CORRECT
%If framemap exceeds number of frames available during tracking (incase of
%behavioural anomily where behav videos were recorded past the experiment)

if length(frameMap1)> length(tracking1(:,1))
    i = 0;
    test = frameMap1(1:find(frameMap1 == length(tracking1(:,1))));
    while isempty(test)
        test = frameMap1(1:find(frameMap1 == length(tracking1(:,1))-i));
        i = i +1;
    end
    frameMap1 = test;
    fprintf('SPECIAL CASE: FrameMap is larger than the behav')
    beep
    pause
end
if length(frameMap2)> length(tracking2(:,1))
    i = 0;
    test = frameMap2(1:find(frameMap2 == length(tracking2(:,1))));
    while isempty(test)
        test = frameMap2(1:find(frameMap2 == length(tracking2(:,1))-i));
        i = i +1;
    end
    frameMap2 = test;
    fprintf('SPECIAL CASE: FrameMap is larger than the behav')
    beep
    pause
end
tracking = cat(1,tracking1(frameMap1,:),tracking2(frameMap2,:));
HD = cat(1,HD1(frameMap1),HD2(frameMap2));

%% Get structure of environment
%Identify where the bounds of the environment are located. Through a subfunction that allows the user to click
%on a plot of the positional data to indicate where the corners of the environment are located.

%% Calculate distances
if isempty(QPOL)
    QPOL = findEdges(tracking);       %finds the edge of the graph and mark the boundaries
end
if isempty(QPOR)
    QPOR = findEdges(tracking);       %finds the edge of the graph and mark the boundaries
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
[disW, ex, ey] = subfuncW(tracking(:,1),tracking(:,2),HD, QPW, degSamp);   %calls funtion to bring back wall distances when neuron fired
dis_rawW = disW;
disW = fillmissing(disW,'pchip',2);                                         %interpolates missing values
disW = disW*pixX;                                                           %Converts boundary distances from pixels to cm.
disW = circshift(disW,90,2);                                                %shifts values by 90 degrees

degSamp = 1;                                                                %angle resolution
[disOL, ex, ey] = subfuncO(tracking(:,1),tracking(:,2),HD, QPOL, degSamp);   %calls funtion to bring back wall distances when neuron fired
dis_rawO = disOL;
disOL = disOL*pixX;                                                           %Converts boundary distances from pixels to cm.
disOL = circshift(disOL,90,2);                                                %shifts values by 90 degrees

degSamp = 1;                                                                %angle resolution
[disOR, ex, ey] = subfuncO(tracking(:,1),tracking(:,2),HD, QPOR, degSamp);   %calls funtion to bring back wall distances when neuron fired
dis_rawO = disOR;
disOR = disOR*pixX;                                                           %Converts boundary distances from pixels to cm.
disOR = circshift(disOR,90,2);                                                %shifts values by 90 degrees
disO = nan(size(disOL));

%Combine Object maps, if conflict take closest object
for i = 1 : length(disOL(:,1))
    for j = 1 : length(disOL(1,:))
        if isnan(disOL(i,j)) || isnan(disOR(i,j))
            disO(i,j)  = nansum([disOL(i,j),disOR(i,j)]);
            if disO(i,j) == 0
                disO(i,j) = NaN;
            end
        elseif disOL(i,j) <= disOR(i,j)
            disO(i,j) = disOL(i,j);
        elseif disOL(i,j) > disOR(i,j)
            disO(i,j) = disOR(i,j);                    
        end        
    end
end

%Loop through every cell, extract and analize firing instances and boundary locations
parfor cellNum = 1 : length(ms.FiltTraces(1,:))
    close all
    mspar = ms;                                                             %create parallel ms variable
    processed = false;                                                      %Progress flag
    stime = tic;                                                            %processing time variable
    distanceBinsPar = distanceBins;                                         %parallel distance variable
    %Exctract calcium trace
    if deconvolve == 1
        fire = mspar.deconvolvedSig(:,cellNum);       
    else
        fire = mspar.FiltTraces(:,cellNum);
    end   
    fire = fire - min(fire);
        
    %Seperate First and second half of recordings (1st vs 2nd recording)
    firsthalflength = length(frameMap1);    
    
    fire1 = fire(1:firsthalflength,1);
    iTime1 = (1:firsthalflength);
    dis1 = disW(1:firsthalflength,:);
    dis1O = disO(1:firsthalflength,:);
    dis1OL = disOL(1:firsthalflength,:);
    dis1OR = disOR(1:firsthalflength,:);
    
    fire2 = fire(firsthalflength+1:end,1);
    iTime2 = (firsthalflength+1:length(tracking(:,1)));
    dis2 = disW((firsthalflength+1:end),:);
    dis2O = disO((firsthalflength+1:end),:);
    dis2OL = disOL((firsthalflength+1:end),:);
    dis2OR = disOR((firsthalflength+1:end),:);
                
    ifire = find(fire);                                                     %Find indices for all non-zero values
    ifire1 = find(fire1);                                                   %Find indices for all non-zero values
    ifire2 = find(fire2);                                                   %Find indices for all non-zero values
    
    if(~isempty(ifire) && ~isempty(ifire1) && ~isempty(ifire2))
        for j = 1 : length(ifire)
            %full run
            ind_fire(j) = ifire(j);                                         %firing index
            cell_x(j) = tracking(((ifire(j))));                     %X postion of mouse during firing
            cell_y(j) = tracking(((ifire(j))),2);                   %Y position of mouse during firing            
            HDfiring(j) = HD((ifire(j)));                           %Head Direction of mouse at time of firing
            %first half
            if j < length(ifire1)
                ind_fire1(j) = ifire1(j);                                   %firing index
                cell_x1(j) = tracking(((iTime1(ifire1(j)))));       %X postion of mouse during Odd min firing
                cell_y1(j) = tracking(((iTime1(ifire1(j)))),2);     %Y position of mouse during Odd min firing                
                HDfiring1(j) = HD((iTime1(ifire1(j))));             %Head Direction of mouse at time of firing
            end
            %second half
            if j < length(ifire2)
                ind_fire2(j) = ifire2(j);                                   %firing index
                cell_x2(j) = tracking(((iTime2(ifire2(j)))));       %X postion of mouse during Even min firing
                cell_y2(j) = tracking(((iTime2(ifire2(j)))),2);     %Y position of mouse during Even min firing               
                HDfiring2(j) = HD((iTime2(ifire2(j))));             %Head Direction of mouse at time of firing
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
        
        occOL = NaN(length(thetaBins), length(distanceBinsPar));                %Object total occupancy bins
        nspkOL = occOL;                                                       %Object occupancy during activity bins
        nspk1OL = nspkOL;
        nspk2OL = nspkOL;
        occ1OL = occOL;
        occ2OL = occOL;
        
        occOR = NaN(length(thetaBins), length(distanceBinsPar));                %Object total occupancy bins
        nspkOR = occOR;                                                       %Object occupancy during activity bins
        nspk1OR = nspkOR;
        nspk2OR = nspkOR;
        occ1OR = occOR;
        occ2OR = occOR;
        
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
            
            tOL = disOL(:,i); %boundary distance for a particular bin
            tO1L = dis1OL(:,i);
            tO2L = dis2OL(:,i);
            
            tOR = disOR(:,i); %boundary distance for a particular bin
            tO1R = dis1OR(:,i);
            tO2R = dis2OR(:,i);
            
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
                
                %Object
                inds = tOL>=distanceBinsPar(k) & tOL<distanceBinsPar(k+1);    %filter through the boundary distances
                occOL(i,k) = sum(inds);                                      %Object occupancy definition
                inds = find(inds);                                          %find all non-zero boundary distances indices
                nspkOL(i,k) = sum(fire(intersect(inds,ci)));                 %Number of spike instances definition
                %first half
                inds1 = tO1L>=distanceBinsPar(k) & tO1L<distanceBinsPar(k+1); %filter through the boundary distances
                occ1OL(i,k) = sum(inds1);                                    %Object occupancy definition
                inds1 = find(inds1);                                        %find all non-zero boundary distances indices
                nspk1OL(i,k) = sum(fire1(intersect(inds1,ci1)));             %Number of spike instances definition
                %second half
                inds2 = tO2L>=distanceBinsPar(k) & tO2L<distanceBinsPar(k+1); %filter through the boundary distances
                occ2OL(i,k) = sum(inds2);                                    %Object occupancy definition
                inds2 = find(inds2);                                        %find all non-zero boundary distances indices
                inds2 = inds2;
                nspk2OL(i,k) = sum(fire2(intersect(inds2,ci2)));             %Number of spike instances definition
                
                %Object
                inds = tOR>=distanceBinsPar(k) & tOR<distanceBinsPar(k+1);    %filter through the boundary distances
                occOR(i,k) = sum(inds);                                      %Object occupancy definition
                inds = find(inds);                                          %find all non-zero boundary distances indices
                nspkOR(i,k) = sum(fire(intersect(inds,ci)));                 %Number of spike instances definition
                %first half
                inds1 = tO1R>=distanceBinsPar(k) & tO1R<distanceBinsPar(k+1); %filter through the boundary distances
                occ1OR(i,k) = sum(inds1);                                    %Object occupancy definition
                inds1 = find(inds1);                                        %find all non-zero boundary distances indices
                nspk1OR(i,k) = sum(fire1(intersect(inds1,ci1)));             %Number of spike instances definition
                %second half
                inds2 = tO2R>=distanceBinsPar(k) & tO2R<distanceBinsPar(k+1); %filter through the boundary distances
                occ2OR(i,k) = sum(inds2);                                    %Object occupancy definition
                inds2 = find(inds2);                                        %find all non-zero boundary distances indices
                inds2 = inds2;
                nspk2OR(i,k) = sum(fire2(intersect(inds2,ci2)));             %Number of spike instances definition
                
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
        
        occ3dOL = zeros(360/DegBinSize,length(occOL(1,:)));
        occ13dOL = zeros(360/DegBinSize,length(occOL(1,:)));
        occ23dOL = zeros(360/DegBinSize,length(occOL(1,:)));
        
        nspk3dOL = zeros(360/DegBinSize,length(occOL(1,:)));
        nspk13dOL = zeros(360/DegBinSize,length(occOL(1,:)));
        nspk23dOL = zeros(360/DegBinSize,length(occOL(1,:)));
        
        occ3dOR = zeros(360/DegBinSize,length(occOR(1,:)));
        occ13dOR = zeros(360/DegBinSize,length(occOR(1,:)));
        occ23dOR = zeros(360/DegBinSize,length(occOR(1,:)));
        
        nspk3dOR = zeros(360/DegBinSize,length(occOR(1,:)));
        nspk13dOR = zeros(360/DegBinSize,length(occOR(1,:)));
        nspk23dOR = zeros(360/DegBinSize,length(occOR(1,:)));
        
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
            
            octempOL = sum(occOL(i*DegBinSize-(DegBinSize-1):i*DegBinSize,:),1);
            octemp1OL = sum(occ1OL(i*DegBinSize-(DegBinSize-1):i*DegBinSize,:),1);
            octemp2OL = sum(occ2OL(i*DegBinSize-(DegBinSize-1):i*DegBinSize,:),1);
            
            occ3dOL(i,:) = octempOL;
            occ13dOL(i,:) = octemp1OL;
            occ23dOL(i,:) = octemp2OL;
            
            nspktempOL = sum(nspkOL(i*DegBinSize-(DegBinSize-1):i*DegBinSize,:),1);
            nspktemp1OL= sum(nspk1OL(i*DegBinSize-(DegBinSize-1):i*DegBinSize,:),1);
            nspktemp2OL = sum(nspk2OL(i*DegBinSize-(DegBinSize-1):i*DegBinSize,:),1);
            
            nspk3dOL(i,:) = nspktempOL;
            nspk13dOL(i,:) = nspktemp1OL;
            nspk23dOL(i,:) = nspktemp2OL;
            
            octempOR = sum(occOR(i*DegBinSize-(DegBinSize-1):i*DegBinSize,:),1);
            octemp1OR = sum(occ1OR(i*DegBinSize-(DegBinSize-1):i*DegBinSize,:),1);
            octemp2OR = sum(occ2OR(i*DegBinSize-(DegBinSize-1):i*DegBinSize,:),1);
            
            occ3dOR(i,:) = octempOR;
            occ13dOR(i,:) = octemp1OR;
            occ23dOR(i,:) = octemp2OR;
            
            nspktempOR = sum(nspkOR(i*DegBinSize-(DegBinSize-1):i*DegBinSize,:),1);
            nspktemp1OR= sum(nspk1OR(i*DegBinSize-(DegBinSize-1):i*DegBinSize,:),1);
            nspktemp2OR = sum(nspk2OR(i*DegBinSize-(DegBinSize-1):i*DegBinSize,:),1);
            
            nspk3dOR(i,:) = nspktempOR;
            nspk13dOR(i,:) = nspktemp1OR;
            nspk23dOR(i,:) = nspktemp2OR;
            
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
        
        occOL = occ3dOL;
        occ1OL = occ13dOL;
        occ2OL = occ23dOL;
        
        nspkOL = nspk3dOL;
        nspk1OL = nspk13dOL;
        nspk2OL = nspk23dOL;
        
        occOR = occ3dOR;
        occ1OR = occ13dOR;
        occ2OR = occ23dOR;
        
        nspkOR = nspk3dOR;
        nspk1OR = nspk13dOR;
        nspk2OR = nspk23dOR;
        
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
            
            %Object
            % bring back to original dims
            occOL=occOL(:,1:end-1);occOL=occOL';
            cutout = find(occOL<0);
            occOL(cutout) = 0;
            nspkOL=nspkOL(:,1:end-1);nspkOL=nspkOL';
            %first half
            occ1OL=occ1OL(:,1:end-1);occ1OL=occ1OL';
            occ1OL(cutout) = 0;
            nspk1OL=nspk1OL(:,1:end-1);nspk1OL=nspk1OL';
            %second half
            occ2OL=occ2OL(:,1:end-1);occ2OL=occ2OL';
            occ2OL(cutout) = 0;
            nspk2OL=nspk2OL(:,1:end-1);nspk2OL=nspk2OL';
            
            rmOL = (nspkOL./occOL) * fps;
            rmOL(find(isnan(rmOL))) = min(rmOL(:));
            rmOL(find(isinf(rmOL))) = min(rmOL(:));
            rmOL = rmOL - min(rmOL(:));
            
            rm1OL = (nspk1OL./occ1OL) * fps;
            rm1OL(find(isnan(rm1OL))) = min(rm1OL(:));
            rm1OL(find(isinf(rm1OL))) = min(rm1OL(:));
            rm1OL = rm1OL - min(rm1OL(:));
            
            rm2OL = (nspk2OL./occ2OL) * fps;
            rm2OL(find(isnan(rm2OL))) = min(rm2OL(:));
            rm2OL(find(isinf(rm2OL))) = min(rm2OL(:));
            rm2OL = rm2OL - min(rm2OL(:));
            
            %Object
            % bring back to original dims
            occOR=occOR(:,1:end-1);occOR=occOR';
            cutout = find(occOR<0);
            occOR(cutout) = 0;
            nspkOR=nspkOR(:,1:end-1);nspkOR=nspkOR';
            %first half
            occ1OR=occ1OR(:,1:end-1);occ1OR=occ1OR';
            occ1OR(cutout) = 0;
            nspk1OR=nspk1OR(:,1:end-1);nspk1OR=nspk1OR';
            %second half
            occ2OR=occ2OR(:,1:end-1);occ2OR=occ2OR';
            occ2OR(cutout) = 0;
            nspk2OR=nspk2OR(:,1:end-1);nspk2OR=nspk2OR';
            
            rmOR = (nspkOR./occOR) * fps;
            rmOR(find(isnan(rmOR))) = min(rmOR(:));
            rmOR(find(isinf(rmOR))) = min(rmOR(:));
            rmOR = rmOR - min(rmOR(:));
            
            rm1OR = (nspk1OR./occ1OR) * fps;
            rm1OR(find(isnan(rm1OR))) = min(rm1OR(:));
            rm1OR(find(isinf(rm1OR))) = min(rm1OR(:));
            rm1OR = rm1OR - min(rm1OR(:));
            
            rm2OR = (nspk2OR./occ2OR) * fps;
            rm2OR(find(isnan(rm2OR))) = min(rm2OR(:));
            rm2OR(find(isinf(rm2OR))) = min(rm2OR(:));
            rm2OR = rm2OR - min(rm2OR(:));
            
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
            corrparO = corr2(rm1OL,rm2O);
            
            %full run
            nd = numel(thetaBins3d);
            rmOL = [rmOL rmOL rmOL];
            rmOL = CMBHOME.Utils.SmoothMat(rmOL,smooth(1:2),smooth(3));   % Smooth it
            rmOL = rmOL(:,nd+1:2*nd); % bring it back
            %first half
            rm1OL = [rm1OL rm1OL rm1OL];
            rm1OL = CMBHOME.Utils.SmoothMat(rm1OL,smooth(1:2),smooth(3));   % Smooth it
            rm1OL = rm1OL(:,nd+1:2*nd); % bring it back
            %second half
            rm2OL = [rm2OL rm2OL rm2OL];
            rm2OL = CMBHOME.Utils.SmoothMat(rm2OL,smooth(1:2),smooth(3));   % Smooth it
            rm2OL = rm2OL(:,nd+1:2*nd); % bring it back
            
            rmOL = fliplr(rmOL);
            rmOL(minDist,:) = 0;
            rm1OL = fliplr(rm1OL);
            rm1OL(minDist,:) = 0;
            rm2OL = fliplr(rm2OL);
            rm2OL(minDist,:) = 0;
            corrparOL = corr2(rm1OL,rm2OL);
            
            %full run
            nd = numel(thetaBins3d);
            rmOR = [rmOR rmOR rmOR];
            rmOR = CMBHOME.Utils.SmoothMat(rmOR,smooth(1:2),smooth(3));   % Smooth it
            rmOR = rmOR(:,nd+1:2*nd); % bring it back
            %first half
            rm1OR = [rm1OR rm1OR rm1OR];
            rm1OR = CMBHOME.Utils.SmoothMat(rm1OR,smooth(1:2),smooth(3));   % Smooth it
            rm1OR = rm1OR(:,nd+1:2*nd); % bring it back
            %second half
            rm2OR = [rm2OR rm2OR rm2OR];
            rm2OR = CMBHOME.Utils.SmoothMat(rm2OR,smooth(1:2),smooth(3));   % Smooth it
            rm2OR = rm2OR(:,nd+1:2*nd); % bring it back
            
            rmOR = fliplr(rmOR);
            rmOR(minDist,:) = 0;
            rm1OR = fliplr(rm1OR);
            rm1OR(minDist,:) = 0;
            rm2OR = fliplr(rm2OR);
            rm2OR(minDist,:) = 0;
            corrparOR = corr2(rm1OR,rm2OR);
            
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
            corrparW = corr2(rm1W,rm2W);
            
            %% EBC METRIC
%             figure;

            %%FIGURE OBJECTS TOGETHER           
            figure('Position', get(0, 'Screensize'));  
%             figure('visible','off');
            n=3;
            m = 3;
            c = 1;
            
            %OBJECT Full run
            % ratemap circular
            subplot(n,m,c);c=c+1;
            % the +pi/2 brings "forwards" to "up"
            [t2, r2] = meshgrid(wrapTo2Pi(thetaBins3d+pi/2), distanceBinsPar(1:end-1));
            [x, y] = pol2cart(t2,r2);
            h=surface(x,y, rmO); shading interp         
            hold on
            set(gca,'XTick',[],'YTick',[])
            axis square
            axis off
            colormap(jet)
            set(gca,'YDir','Normal')
            freezeColors
            title('rm both Obj')
            info = ['Max: ' num2str(max(max(rmO)))];
            text(30,8,info);            
            
            %Trajectory map
            if deconvolve == 1
                subplot(n,m,c);c=c+1;
%                 edg = splitter(QPO);
                hold on
                plot(pixX*tracking(:,1),-pixY*tracking(:,2),'Color',[.7 .7 .7])
                colormap(hsv)
                xlim(pixX*[min(tracking(:,1)) max(tracking(:,1))]);ylim(pixY*[-max(tracking(:,2)) -min(tracking(:,2))])
                cx=pixX*cell_x(:);
                cy=pixY*cell_y(:);
                axis off
                axis square
                scatter(cx,-cy,38,HDfiring(:),'filled')
                set(gca,'YDir','Normal')
                caxis([0 360])
                title('Traj')                
            else
                subplot(n,m,c);c=c+1;
%                 edg = splitter(QPO);
                hold on
                plot(pixX*tracking(:,1),-pixY*tracking(:,2),'Color',[.7 .7 .7])
                colormap(hsv) 
                xlim(pixX*[min(tracking(:,1)) max(tracking(:,1))]);ylim(pixY*[-max(tracking(:,2)) -min(tracking(:,2))])                               
            end
            set(gca,'YDir','Normal')    
            caxis([0 360])
            title('Traj')
            axis off
            axis square
            
            %Object EBC Metric             
            metricO = mean(rmO,1)';
            if ~deconvolve
                metricO = metricO - min(metricO);
            end
            
            xs = metricO(1:end-1).*cos(degBins(1:end-1)); % average
            ys = metricO(1:end-1).*sin(degBins(1:end-1));                       
            
            ang_hdO = atan2(mean(ys),mean(xs)); % mean direction
            
            mrO = (cos(ang_hdO)*sum(xs) + sin(ang_hdO)*sum(ys)) / sum(metricO(1:end-1)); % mean resultant length
            
            %Object polar plot      
            subplot(n,m,c);c=c+1;
            polarplot(degBins,metricO)
            coordlims=axis;
            hold on;
            polarplot([ang_hdO ang_hdO ],[0 mrO], 'r')
            pol = gca;
            pol.ThetaZeroLocation = 'top';
            pol.ThetaColor = 'none';
            pol.RTickLabel = [];
            hold off
            title('Object Directionality')
            stat = ['MRL: ' num2str(mrO) newline 'Angle : ' num2str(rad2deg(ang_hdO))];
            text(0.2,coordlims(4),stat);                       
            
            %FIRST HALF
            % Object ratemap circular
            subplot(n,m,c);c=c+1;
            % the +pi/2 brings "forwards" to "up"
            [t2, r2] = meshgrid(wrapTo2Pi(thetaBins3d+pi/2), distanceBinsPar(1:end-1));
            [x, y] = pol2cart(t2,r2);
            h=surface(x,y, rm1O); shading interp
            hold on
            set(gca,'XTick',[],'YTick',[])
            axis square
            axis off
            colormap(jet)
            set(gca,'YDir','Normal')
            freezeColors
            title('Object 1st Half Ratemap')
            info = ['Max: ' num2str(max(max(rm1O))) newline ' Correlation: ' num2str(corrparO)];
            text(30,8,info);                        

            %Trajectory map
            subplot(n,m,c);c=c+1;
            if deconvolve == 1
%                 edg = splitter(QP);
                hold on
                plot(pixX*tracking((iTime1),1),-pixY*tracking((iTime1),2),'Color',[.7 .7 .7])
                colormap(hsv)
                xlim(pixX*[min(tracking(:,1)) max(tracking(:,1))]);ylim(pixY*[-max(tracking(:,2)) -min(tracking(:,2))])
                cx=pixX*cell_x1(:);
                cy=pixY*cell_y1(:);
                axis off
                axis square
                scatter(cx,-cy,38,HDfiring1(:),'filled')
                set(gca,'YDir','Normal')
                caxis([0 360])
                title('1st Half Trajectory')                
            else
%                 edg = splitter(QP);
                hold on
                plot(pixX*tracking((iTime1),1),-pixY*tracking((iTime1),2),'Color',[.7 .7 .7])
                colormap(hsv)
                xlim(pixX*[min(tracking(:,1)) max(tracking(:,1))]);ylim(pixY*[-max(tracking(:,2)) -min(tracking(:,2))])
                set(gca,'YDir','Normal')
                title('1st Half Trajectory')                 
            end
            set(gca,'YDir','Normal')    
            caxis([0 360])
            title('1st Half Trajectory')
            axis off
            axis square                        
            
            %Object EBC Metric
            metric1O = mean(rm1O,1)';
            if ~deconvolve
                metric1O = metric1O - min(metric1O);
            end
            
            xs = metric1O(1:end-1).*cos(degBins(1:end-1)); % average
            ys = metric1O(1:end-1).*sin(degBins(1:end-1));                        
            
            ang_hd1O = atan2(mean(ys),mean(xs)); % mean direction
            
            mr1O = (cos(ang_hd1O)*sum(xs) + sin(ang_hd1O)*sum(ys)) / sum(metric1O(1:end-1)); % mean resultant length
            
            %Object Polar plot
            subplot(n,m,c);c=c+1;
            polarplot(degBins,metric1O)
            coordlims=axis;
            hold on;
            polarplot([ang_hd1O ang_hd1O ],[0 mr1O], 'r')
            pol = gca;
            pol.ThetaZeroLocation = 'top';
            pol.ThetaColor = 'none';
            pol.RTickLabel = [];
            hold off
            title('Object Directionality: 1st Half')
            stat = ['MRL: ' num2str(mr1O) newline 'Angle : ' num2str(rad2deg(ang_hd1O))];
            text(0.2,coordlims(4),stat);                        
            
            %SECOND HALF
            % Object ratemap circular
            subplot(n,m,c);c=c+1;
            % the +pi/2 brings "forwards" to "up"
            [t2, r2] = meshgrid(wrapTo2Pi(thetaBins3d+pi/2), distanceBinsPar(1:end-1));
            [x, y] = pol2cart(t2,r2);
            h=surface(x,y, rm2O); shading interp
            hold on
            set(gca,'XTick',[],'YTick',[])
            axis square
            axis off
            colormap(jet)
            set(gca,'YDir','Normal')
            freezeColors
            title('Object 2nd Half Ratemap')
            info = ['Max: ' num2str(max(max(rm2O)))];
            text(30,8,info);                        

            %Trajectory map
            subplot(n,m,c);c=c+1;
            if deconvolve == 1
%                 edg = splitter(QP);
                hold on
                plot(pixX*tracking((iTime2),1),-pixY*tracking((iTime2),2),'Color',[.7 .7 .7])
                colormap(hsv)
                xlim(pixX*[min(tracking(:,1)) max(tracking(:,1))]);ylim(pixY*[-max(tracking(:,2)) -min(tracking(:,2))])
                cx=pixX*cell_x2(:);
                cy=pixY*cell_y2(:);
                axis off
                axis square
                scatter(cx,-cy,38,HDfiring2(:),'filled')
                set(gca,'YDir','Normal')
                caxis([0 360])
                title('2nd Half Trajectory')                
            else
%                 edg = splitter(QP);
                hold on
                plot(pixX*tracking((iTime2),1),-pixY*tracking((iTime2),2),'Color',[.7 .7 .7])
                colormap(hsv)
                xlim(pixX*[min(tracking(:,1)) max(tracking(:,1))]);ylim(pixY*[-max(tracking(:,2)) -min(tracking(:,2))])
                set(gca,'YDir','Normal')
                title('2nd Half Trajectory')                 
            end
            set(gca,'YDir','Normal')    
            caxis([0 360])
            title('2nd Half Trajectory')
            axis off
            axis square                        
            
            %Object EBC Metric
            metric2O = mean(rm2O,1)';
            if ~deconvolve
                metric2O = metric2O - min(metric2O);
            end
            
            xs = metric2O(1:end-1).*cos(degBins(1:end-1)); % average
            ys = metric2O(1:end-1).*sin(degBins(1:end-1));
            
            ang_hd2O = atan2(mean(ys),mean(xs)); % mean direction
            
            mr2O = (cos(ang_hd2O)*sum(xs) + sin(ang_hd2O)*sum(ys)) / sum(metric2O(1:end-1)); % mean resultant length
            
            %Object Polar plot
            subplot(n,m,c);c=c+1;
            polarplot(degBins,metric2O)
            coordlims=axis;
            hold on;
            polarplot([ang_hd2O ang_hd2O ],[0 mr2O], 'r')
            pol = gca;
            pol.ThetaZeroLocation = 'top';
            pol.ThetaColor = 'none';
            pol.RTickLabel = [];
            hold off
            title('Ojbect Directionality: 2nd Half')
            stat = ['MRL: ' num2str(mr2O) newline 'Angle : ' num2str(rad2deg(ang_hd2O))];
            text(0.2,coordlims(4),stat);
            
            saveas(gcf,[pwd '/' name,'/ObjectsTogether/Objects',num2str(cellNum),'EBC.jpg']);
            
            
            %FIGURE OBJECTS INDIVIDUAL
            figure('Position', get(0, 'Screensize'));  
%             figure('visible','off');
            n=3;
            m = 5;
            c = 1;
            
            %OBJECT1 Full run
            % ratemap circular
            subplot(n,m,c);c=c+1;
            % the +pi/2 brings "forwards" to "up"
            [t2, r2] = meshgrid(wrapTo2Pi(thetaBins3d+pi/2), distanceBinsPar(1:end-1));
            [x, y] = pol2cart(t2,r2);
            h=surface(x,y, rmOL); shading interp         
            hold on
            set(gca,'XTick',[],'YTick',[])
            axis square
            axis off
            colormap(jet)
            set(gca,'YDir','Normal')
            freezeColors
            title('rm Left Obj')
            info = ['Max: ' num2str(max(max(rmOL)))];
            text(30,8,info);
            
            %OBJECT2 Full run
            % ratemap circular
            subplot(n,m,c);c=c+1;
            % the +pi/2 brings "forwards" to "up"
            [t2, r2] = meshgrid(wrapTo2Pi(thetaBins3d+pi/2), distanceBinsPar(1:end-1));
            [x, y] = pol2cart(t2,r2);
            h=surface(x,y, rmOR); shading interp         
            hold on
            set(gca,'XTick',[],'YTick',[])
            axis square
            axis off
            colormap(jet)
            set(gca,'YDir','Normal')
            freezeColors
            title('rm Right Obj')
            info = ['Max: ' num2str(max(max(rmOR)))];
            text(30,8,info);
            
            %Trajectory map
            if deconvolve == 1
                subplot(n,m,c);c=c+1;
%                 edg = splitter(QPO);
                hold on
                plot(pixX*tracking(:,1),-pixY*tracking(:,2),'Color',[.7 .7 .7])
                colormap(hsv)
                xlim(pixX*[min(tracking(:,1)) max(tracking(:,1))]);ylim(pixY*[-max(tracking(:,2)) -min(tracking(:,2))])
                cx=pixX*cell_x(:);
                cy=pixY*cell_y(:);
                axis off
                axis square
                scatter(cx,-cy,38,HDfiring(:),'filled')
                set(gca,'YDir','Normal')
                caxis([0 360])
                title('Traj')                
            else
                subplot(n,m,c);c=c+1;
%                 edg = splitter(QPO);
                hold on
                plot(pixX*tracking(:,1),-pixY*tracking(:,2),'Color',[.7 .7 .7])
                colormap(hsv) 
                xlim(pixX*[min(tracking(:,1)) max(tracking(:,1))]);ylim(pixY*[-max(tracking(:,2)) -min(tracking(:,2))])                               
            end
            set(gca,'YDir','Normal')    
            caxis([0 360])
            title('Traj')
            axis off
            axis square
            
            %Object EBC Metric             
            metricOL = mean(rmOL,1)';
            if ~deconvolve
                metricOL = metricOL - min(metricOL);
            end
            
            xs = metricOL(1:end-1).*cos(degBins(1:end-1)); % average
            ys = metricOL(1:end-1).*sin(degBins(1:end-1));                       
            
            ang_hdOL = atan2(mean(ys),mean(xs)); % mean direction
            
            mrOL = (cos(ang_hdOL)*sum(xs) + sin(ang_hdOL)*sum(ys)) / sum(metricOL(1:end-1)); % mean resultant length
            
            %Object polar plot      
            subplot(n,m,c);c=c+1;
            polarplot(degBins,metricOL)
            coordlims=axis;
            hold on;
            polarplot([ang_hdOL ang_hdOL ],[0 mrOL], 'r')
            pol = gca;
            pol.ThetaZeroLocation = 'top';
            pol.ThetaColor = 'none';
            pol.RTickLabel = [];
            hold off
            title('Left Object Directionality')
            stat = ['MRL: ' num2str(mrOL) newline 'Angle : ' num2str(rad2deg(ang_hdOL))];
            text(0.2,coordlims(4),stat);
            
            %Object 2 EBC Metric             
            metricOR = mean(rmOR,1)';
            if ~deconvolve
                metricOR = metricOR - min(metricOR);
            end
            
            xs = metricOR(1:end-1).*cos(degBins(1:end-1)); % average
            ys = metricOR(1:end-1).*sin(degBins(1:end-1));                       
            
            ang_hdOR = atan2(mean(ys),mean(xs)); % mean direction
            
            mrOR = (cos(ang_hdOR)*sum(xs) + sin(ang_hdOR)*sum(ys)) / sum(metricOR(1:end-1)); % mean resultant length
            
            %Object polar plot      
            subplot(n,m,c);c=c+1;
            polarplot(degBins,metricOR)
            coordlims=axis;
            hold on;
            polarplot([ang_hdOR ang_hdOR ],[0 mrOR], 'r')
            pol = gca;
            pol.ThetaZeroLocation = 'top';
            pol.ThetaColor = 'none';
            pol.RTickLabel = [];
            hold off
            title('Right Object Directionality')
            stat = ['MRL: ' num2str(mrOR) newline 'Angle : ' num2str(rad2deg(ang_hdOR))];
            text(0.2,coordlims(4),stat);
            
            %FIRST HALF
            % Object 1 ratemap circular
            subplot(n,m,c);c=c+1;
            % the +pi/2 brings "forwards" to "up"
            [t2, r2] = meshgrid(wrapTo2Pi(thetaBins3d+pi/2), distanceBinsPar(1:end-1));
            [x, y] = pol2cart(t2,r2);
            h=surface(x,y, rm1OL); shading interp
            hold on
            set(gca,'XTick',[],'YTick',[])
            axis square
            axis off
            colormap(jet)
            set(gca,'YDir','Normal')
            freezeColors
            title('Object 1st half Ratemap')
            info = ['Max: ' num2str(max(max(rm1OL))) newline ' Correlation: ' num2str(corrparOL)];
            text(30,8,info);
                        
            % Object 2 ratemap circular
            subplot(n,m,c);c=c+1;
            % the +pi/2 brings "forwards" to "up"
            [t2, r2] = meshgrid(wrapTo2Pi(thetaBins3d+pi/2), distanceBinsPar(1:end-1));
            [x, y] = pol2cart(t2,r2);
            h=surface(x,y, rm1OR); shading interp
            hold on
            set(gca,'XTick',[],'YTick',[])
            axis square
            axis off
            colormap(jet)
            set(gca,'YDir','Normal')
            freezeColors
            title('Object 1st Half Ratemap')
            info = ['Max: ' num2str(max(max(rm1OR))) newline ' Correlation: ' num2str(corrparOR)];
            text(30,8,info);

            %Trajectory map
            subplot(n,m,c);c=c+1;
            if deconvolve == 1
%                 edg = splitter(QP);
                hold on
                plot(pixX*tracking((iTime1),1),-pixY*tracking((iTime1),2),'Color',[.7 .7 .7])
                colormap(hsv)
                xlim(pixX*[min(tracking(:,1)) max(tracking(:,1))]);ylim(pixY*[-max(tracking(:,2)) -min(tracking(:,2))])
                cx=pixX*cell_x1(:);
                cy=pixY*cell_y1(:);
                axis off
                axis square
                scatter(cx,-cy,38,HDfiring1(:),'filled')
                set(gca,'YDir','Normal')
                caxis([0 360])
                title('1st Half Trajectory')                
            else
%                 edg = splitter(QP);
                hold on
                plot(pixX*tracking((iTime1),1),-pixY*tracking((iTime1),2),'Color',[.7 .7 .7])
                colormap(hsv)
                xlim(pixX*[min(tracking(:,1)) max(tracking(:,1))]);ylim(pixY*[-max(tracking(:,2)) -min(tracking(:,2))])
                set(gca,'YDir','Normal')
                title('1st Half Trajectory')                 
            end
            set(gca,'YDir','Normal')    
            caxis([0 360])
            title('1st Half Trajectory')
            axis off
            axis square                        
            
            %Object 1 EBC Metric
            metric1OL = mean(rm1OL,1)';
            if ~deconvolve
                metric1OL = metric1OL - min(metric1OL);
            end
            
            xs = metric1OL(1:end-1).*cos(degBins(1:end-1)); % average
            ys = metric1OL(1:end-1).*sin(degBins(1:end-1));                        
            
            ang_hd1OL = atan2(mean(ys),mean(xs)); % mean direction
            
            mr1OL = (cos(ang_hd1OL)*sum(xs) + sin(ang_hd1OL)*sum(ys)) / sum(metric1OL(1:end-1)); % mean resultant length
            
            %Object Polar plot
            subplot(n,m,c);c=c+1;
            polarplot(degBins,metric1OL)
            coordlims=axis;
            hold on;
            polarplot([ang_hd1OL ang_hd1OL ],[0 mr1OL], 'r')
            pol = gca;
            pol.ThetaZeroLocation = 'top';
            pol.ThetaColor = 'none';
            pol.RTickLabel = [];
            hold off
            title('Left Object Directionality: 1st Half')
            stat = ['MRL: ' num2str(mr1OL) newline 'Angle : ' num2str(rad2deg(ang_hd1OL))];
            text(0.2,coordlims(4),stat);
            
            %Object 2 EBC Metric
            metric1OR = mean(rm1OR,1)';
            if ~deconvolve
                metric1OR = metric1OR - min(metric1OR);
            end
            
            xs = metric1OR(1:end-1).*cos(degBins(1:end-1)); % average
            ys = metric1OR(1:end-1).*sin(degBins(1:end-1));                        
            
            ang_hd1OR = atan2(mean(ys),mean(xs)); % mean direction
            
            mr1OR = (cos(ang_hd1OR)*sum(xs) + sin(ang_hd1OR)*sum(ys)) / sum(metric1OR(1:end-1)); % mean resultant length
            
            %Object Polar plot
            subplot(n,m,c);c=c+1;
            polarplot(degBins,metric1OR)
            coordlims=axis;
            hold on;
            polarplot([ang_hd1OR ang_hd1OR],[0 mr1OR], 'r')
            pol = gca;
            pol.ThetaZeroLocation = 'top';
            pol.ThetaColor = 'none';
            pol.RTickLabel = [];
            hold off
            title('Right Object Directionality: 1st Half')
            stat = ['MRL: ' num2str(mr1OR) newline 'Angle : ' num2str(rad2deg(ang_hd1OR))];
            text(0.2,coordlims(4),stat);
            
            %SECOND HALF
            % Object 1 ratemap circular
            subplot(n,m,c);c=c+1;
            % the +pi/2 brings "forwards" to "up"
            [t2, r2] = meshgrid(wrapTo2Pi(thetaBins3d+pi/2), distanceBinsPar(1:end-1));
            [x, y] = pol2cart(t2,r2);
            h=surface(x,y, rm2OL); shading interp
            hold on
            set(gca,'XTick',[],'YTick',[])
            axis square
            axis off
            colormap(jet)
            set(gca,'YDir','Normal')
            freezeColors
            title('Left Object 2nd Half Ratemap')
            info = ['Max: ' num2str(max(max(rm2OL)))];
            text(30,8,info);
            
            % Object ratemap circular
            subplot(n,m,c);c=c+1;
            % the +pi/2 brings "forwards" to "up"
            [t2, r2] = meshgrid(wrapTo2Pi(thetaBins3d+pi/2), distanceBinsPar(1:end-1));
            [x, y] = pol2cart(t2,r2);
            h=surface(x,y, rm2OR); shading interp
            hold on
            set(gca,'XTick',[],'YTick',[])
            axis square
            axis off
            colormap(jet)
            set(gca,'YDir','Normal')
            freezeColors
            title('Right Object 2nd Half Ratemap')
            info = ['Max: ' num2str(max(max(rm2OR)))];
            text(30,8,info);

            %Trajectory map
            subplot(n,m,c);c=c+1;
            if deconvolve == 1
%                 edg = splitter(QP);
                hold on
                plot(pixX*tracking(iTime2,1),-pixY*tracking(iTime2,2),'Color',[.7 .7 .7])
                colormap(hsv)
                xlim(pixX*[min(tracking(:,1)) max(tracking(:,1))]);ylim(pixY*[-max(tracking(:,2)) -min(tracking(:,2))])
                cx=pixX*cell_x2(:);
                cy=pixY*cell_y2(:);
                axis off
                axis square
                scatter(cx,-cy,38,HDfiring2(:),'filled')
                set(gca,'YDir','Normal')
                caxis([0 360])
                title('2nd Half Trajectory')                
            else
%                 edg = splitter(QP);
                hold on
                plot(pixX*tracking((iTime2),1),-pixY*tracking((iTime2),2),'Color',[.7 .7 .7])
                colormap(hsv)
                xlim(pixX*[min(tracking(:,1)) max(tracking(:,1))]);ylim(pixY*[-max(tracking(:,2)) -min(tracking(:,2))])
                set(gca,'YDir','Normal')
                title('2nd Half Trajectory')                 
            end
            set(gca,'YDir','Normal')    
            caxis([0 360])
            title('2nd Half Trajectory')
            axis off
            axis square                        
            
            %Object 1 EBC Metric
            metric2OL = mean(rm2OL,1)';
            if ~deconvolve
                metric2OL = metric2OL - min(metric2OL);
            end
            
            xs = metric2OL(1:end-1).*cos(degBins(1:end-1)); % average
            ys = metric2OL(1:end-1).*sin(degBins(1:end-1));
            
            ang_hd2OL = atan2(mean(ys),mean(xs)); % mean direction
            
            mr2OL = (cos(ang_hd2OL)*sum(xs) + sin(ang_hd2OL)*sum(ys)) / sum(metric2OL(1:end-1)); % mean resultant length
            
            %Object Polar plot
            subplot(n,m,c);c=c+1;
            polarplot(degBins,metric2OL)
            coordlims=axis;
            hold on;
            polarplot([ang_hd2OL ang_hd2OL ],[0 mr2OL], 'r')
            pol = gca;
            pol.ThetaZeroLocation = 'top';
            pol.ThetaColor = 'none';
            pol.RTickLabel = [];
            hold off
            title('Left Ojbect Directionality: 2nd Half')
            stat = ['MRL: ' num2str(mr2OL) newline 'Angle : ' num2str(rad2deg(ang_hd2OL))];
            text(0.2,coordlims(4),stat);
            
            %Object EBC Metric
            metric2OR = mean(rm2OR,1)';
            if ~deconvolve
                metric2OR = metric2OR - min(metric2OR);
            end
            
            xs = metric2OR(1:end-1).*cos(degBins(1:end-1)); % average
            ys = metric2OR(1:end-1).*sin(degBins(1:end-1));
            
            ang_hd2OR = atan2(mean(ys),mean(xs)); % mean direction
            
            mr2OR = (cos(ang_hd2OR)*sum(xs) + sin(ang_hd2OR)*sum(ys)) / sum(metric2OR(1:end-1)); % mean resultant length
            
            %Object Polar plot
            subplot(n,m,c);c=c+1;
            polarplot(degBins,metric2OR)
            coordlims=axis;
            hold on;
            polarplot([ang_hd2OR ang_hd2OR ],[0 mr2OR], 'r')
            pol = gca;
            pol.ThetaZeroLocation = 'top';
            pol.ThetaColor = 'none';
            pol.RTickLabel = [];
            hold off
            title('Right Ojbect Directionality: 2nd Half')
            stat = ['MRL: ' num2str(mr2OR) newline 'Angle : ' num2str(rad2deg(ang_hd2OR))];
            text(0.2,coordlims(4),stat);
            
            saveas(gcf,[name,'/ObjectSeperate/Seperate',num2str(cellNum),'EBC.jpg']);
            
            %FIGURE WALLS
            figure('Position', get(0, 'Screensize'));  
%             figure('visible','off');
            n=3;
            m = 3;
            c = 1;
            
            %Wall Full run
            % ratemap circular
            subplot(n,m,c);c=c+1;
            % the +pi/2 brings "forwards" to "up"
            [t2, r2] = meshgrid(wrapTo2Pi(thetaBins3d+pi/2), distanceBinsPar(1:end-1));
            [x, y] = pol2cart(t2,r2);
            h=surface(x,y, rmW); shading interp                   
            hold on
            set(gca,'XTick',[],'YTick',[])
            axis square
            axis off
            colormap(jet)
            set(gca,'YDir','Normal')
            freezeColors
            title('rm Wall')
            info = ['Max: ' num2str(max(max(rmW)))];
            text(30,8,info);
            
            %Trajectory map
            if deconvolve == 1
                subplot(n,m,c);c=c+1;
%                 edg = splitter(QPO);
                hold on
                plot(pixX*tracking(:,1),-pixY*tracking(:,2),'Color',[.7 .7 .7])
                colormap(hsv)
                xlim(pixX*[min(tracking(:,1)) max(tracking(:,1))]);ylim(pixY*[-max(tracking(:,2)) -min(tracking(:,2))])
                cx=pixX*cell_x(:);
                cy=pixY*cell_y(:);
                axis off
                axis square
                scatter(cx,-cy,38,HDfiring(:),'filled')
                set(gca,'YDir','Normal')
                caxis([0 360])
                title('Traj')                
            else
                subplot(n,m,c);c=c+1;
%                 edg = splitter(QPO);
                hold on
                plot(pixX*tracking(:,1),-pixY*tracking(:,2),'Color',[.7 .7 .7])
                colormap(hsv) 
                xlim(pixX*[min(tracking(:,1)) max(tracking(:,1))]);ylim(pixY*[-max(tracking(:,2)) -min(tracking(:,2))])                               
            end
            set(gca,'YDir','Normal')    
            caxis([0 360])
            title('Traj')
            axis off
            axis square                       
            
            %Wall EBC Metric
            metricW = mean(rmW,1)';
            if ~deconvolve
                metricW = metricW - min(metricW);
            end
            
            xs = metricW(1:end-1).*cos(degBins(1:end-1)); % average
            ys = metricW(1:end-1).*sin(degBins(1:end-1));
            
            ang_hdW = atan2(mean(ys),mean(xs)); % mean direction
            
            mrW = (cos(ang_hdW)*sum(xs) + sin(ang_hdW)*sum(ys)) / sum(metricW(1:end-1)); % mean resultant length
            
            %Wall polar plot
            subplot(n,m,c);c=c+1;
            polarplot(degBins,metricW)
            coordlims=axis;
            hold on;
            polarplot([ang_hdW ang_hdW ],[0 mrW], 'r')
            pol = gca;
            pol.ThetaZeroLocation = 'top';
            pol.ThetaColor = 'none';
            pol.RTickLabel = [];
            hold off
            title('Wall Directionality')
            stat = ['MRL: ' num2str(mrW) newline 'Angle : ' num2str(rad2deg(ang_hdW))];
            text(0.2,coordlims(4),stat);
            
            %FIRST HALF            
            
            % Wall ratemap circular
            subplot(n,m,c);c=c+1;
            % the +pi/2 brings "forwards" to "up"
            [t2, r2] = meshgrid(wrapTo2Pi(thetaBins3d+pi/2), distanceBinsPar(1:end-1));
            [x, y] = pol2cart(t2,r2);
            h=surface(x,y, rm1W); shading interp
            hold on
            set(gca,'XTick',[],'YTick',[])
            axis square
            axis off
            colormap(jet)
            set(gca,'YDir','Normal')
            freezeColors
            title('Wall 1st Half Ratemap')
            info = ['Max: ' num2str(max(max(rm1W))) newline ' Correlation: ' num2str(corrparW)];
            text(30,8,info);

            %Trajectory map
            subplot(n,m,c);c=c+1;
            if deconvolve == 1
%                 edg = splitter(QP);
                hold on
                plot(pixX*tracking((iTime1),1),-pixY*tracking((iTime1),2),'Color',[.7 .7 .7])
                colormap(hsv)
                xlim(pixX*[min(tracking(:,1)) max(tracking(:,1))]);ylim(pixY*[-max(tracking(:,2)) -min(tracking(:,2))])
                cx=pixX*cell_x1(:);
                cy=pixY*cell_y1(:);
                axis off
                axis square
                scatter(cx,-cy,38,HDfiring1(:),'filled')
                set(gca,'YDir','Normal')
                caxis([0 360])
                title('1st Half Trajectory')                
            else
%                 edg = splitter(QP);
                hold on
                plot(pixX*tracking((iTime1),1),-pixY*tracking((iTime1),2),'Color',[.7 .7 .7])
                colormap(hsv)
                xlim(pixX*[min(tracking(:,1)) max(tracking(:,1))]);ylim(pixY*[-max(tracking(:,2)) -min(tracking(:,2))])
                set(gca,'YDir','Normal')
                title('1st Half Trajectory')                 
            end
            set(gca,'YDir','Normal')    
            caxis([0 360])
            title('1st Half Trajectory')
            axis off
            axis square                                               
            
            %Wall EBC Metric
            metric1W = mean(rm1W,1)';
            if ~deconvolve
                metric1W = metric1W - min(metric1W);
            end
            
            xs = metric1W(1:end-1).*cos(degBins(1:end-1)); % average
            ys = metric1W(1:end-1).*sin(degBins(1:end-1));                        
            
            ang_hd1W = atan2(mean(ys),mean(xs)); % mean direction
            
            mr1W = (cos(ang_hd1W)*sum(xs) + sin(ang_hd1W)*sum(ys)) / sum(metric1W(1:end-1)); % mean resultant length
            
            %Wall Polar plot
            subplot(n,m,c);c=c+1;
            polarplot(degBins,metric1W)
            coordlims=axis;
            hold on;
            polarplot([ang_hd1W ang_hd1W ],[0 mr1W], 'r')
            pol = gca;
            pol.ThetaZeroLocation = 'top';
            pol.ThetaColor = 'none';
            pol.RTickLabel = [];
            hold off
            title('Wall Directionality: 1st Half')
            stat = ['MRL: ' num2str(mr1W) newline 'Angle : ' num2str(rad2deg(ang_hd1W))];
            text(0.2,coordlims(4),stat);
            
            %SECOND HALF            
            % Wall ratemap circular
            subplot(n,m,c);c=c+1;
            % the +pi/2 brings "forwards" to "up"
            [t2, r2] = meshgrid(wrapTo2Pi(thetaBins3d+pi/2), distanceBinsPar(1:end-1));
            [x, y] = pol2cart(t2,r2);
            h=surface(x,y, rm2W); shading interp
            hold on
            set(gca,'XTick',[],'YTick',[])
            axis square
            axis off
            colormap(jet)
            set(gca,'YDir','Normal')
            freezeColors
            title('Wall 2nd Half Ratemap')
            info = ['Max: ' num2str(max(max(rm2W)))];
            text(30,8,info);

            %Trajectory map
            subplot(n,m,c);c=c+1;
            if deconvolve == 1
%                 edg = splitter(QP);
                hold on
                plot(pixX*tracking((iTime2),1),-pixY*tracking((iTime2),2),'Color',[.7 .7 .7])
                colormap(hsv)
                xlim(pixX*[min(tracking(:,1)) max(tracking(:,1))]);ylim(pixY*[-max(tracking(:,2)) -min(tracking(:,2))])
                cx=pixX*cell_x2(:);
                cy=pixY*cell_y2(:);
                axis off
                axis square
                scatter(cx,-cy,38,HDfiring2(:),'filled')
                set(gca,'YDir','Normal')
                caxis([0 360])
                title('2nd Half Trajectory')                
            else
%                 edg = splitter(QP);
                hold on
                plot(pixX*tracking((iTime2),1),-pixY*tracking((iTime2),2),'Color',[.7 .7 .7])
                colormap(hsv)
                xlim(pixX*[min(tracking(:,1)) max(tracking(:,1))]);ylim(pixY*[-max(tracking(:,2)) -min(tracking(:,2))])
                set(gca,'YDir','Normal')
                title('2nd Half Trajectory')                 
            end
            set(gca,'YDir','Normal')    
            caxis([0 360])
            title('2nd Half Trajectory')
            axis off
            axis square                        
                                   
            %Wall EBC Metric
            metric2W = mean(rm2W,1)';
            if ~deconvolve
                metric2W = metric2W - min(metric2W);
            end
            
            xs = metric2W(1:end-1).*cos(degBins(1:end-1)); % average
            ys = metric2W(1:end-1).*sin(degBins(1:end-1));                        
            
            ang_hd2W = atan2(mean(ys),mean(xs)); % mean direction
            
            mr2W = (cos(ang_hd2W)*sum(xs) + sin(ang_hd2W)*sum(ys)) / sum(metric2W(1:end-1)); % mean resultant length
            
            %Wall Polar plot
            subplot(n,m,c);c=c+1;
            polarplot(degBins,metric1W)
            coordlims=axis;
            hold on;
            polarplot([ang_hd2W ang_hd2W ],[0 mr2W], 'r')
            pol = gca;
            pol.ThetaZeroLocation = 'top';
            pol.ThetaColor = 'none';
            pol.RTickLabel = [];
            hold off
            title('Wall Directionality: 2nd Half')
            stat = ['MRL: ' num2str(mr2W) newline 'Angle : ' num2str(rad2deg(ang_hd2W))];
            text(0.2,coordlims(4),stat);
            
            %Save Results
            ind_fire = []; %Indices of neuron activity/firing
            ind_fire1 = [];
            ind_fire2 =[];
            
            mrlittO(cellNum) = mrO; %save the MRL
            mrlittO1(cellNum) =mr1O;
            mrlittO2(cellNum) =mr2O;
            corrittO(cellNum) = corrparO;
            
            mrlittOL(cellNum) = mrOL; %save the MRL
            mrlittO1L(cellNum) =mr1OL;
            mrlittO2L(cellNum) =mr2OL;
            corrittOL(cellNum) = corrparOL;
            
            mrlittOR(cellNum) = mrOR; %save the MRL
            mrlittO1R(cellNum) =mr1OR;
            mrlittO2R(cellNum) =mr2OR;
            corrittOR(cellNum) = corrparOR;
            
            mrlittW(cellNum) = mrW; %save the MRL
            mrlittW1(cellNum) =mr1W;
            mrlittW2(cellNum) =mr2W;
            corrittW(cellNum) = corrparW;     
            
            saveas(gcf,[name,'/wall/wall',num2str(cellNum),'EBC.jpg']); %saving figure as a picture file (.jpg) in the new folder "EBCresults"
            
            clf
            ind_fire = NaN(ms.numFrames,length(ms.FiltTraces(1,:)));        %Indices of neuron activity/firing
            cell_x = NaN(ms.numFrames,length(ms.FiltTraces(1,:)));          %X coordinate of locations where cell fired
            cell_y = NaN(ms.numFrames,length(ms.FiltTraces(1,:)));          %Y cooridnates of locations where cell fired
            HDfiring = ind_fire;                                            %Indices of activity for head direction 
            ind_fire1 = NaN(size(ind_fire1));
            cell_x1 = NaN;
            cell_y1 = NaN;
            HDfiring1 = NaN;
            ind_fire2 = NaN(size(ind_fire2));
            cell_x2 = NaN;
            cell_y2 = NaN;
            HDfiring2 = NaN;
            
            processed = true;
        end
    end
    etime = toc(stime);
    timevar(cellNum) = etime;
    if processed
        ratemapsO(:,:,cellNum) = rmO;        
        ratemaps1O(:,:,cellNum) = rm1O;        
        ratemaps2O(:,:,cellNum) = rm2O;
        headangleO(1,cellNum) = ang_hdO;
        headangle1O(1,cellNum) = ang_hd1O;
        headangle2O(1,cellNum) = ang_hd2O;
        
        ratemapsOL(:,:,cellNum) = rmOL;        
        ratemaps1OL(:,:,cellNum) = rm1OL;        
        ratemaps2OL(:,:,cellNum) = rm2OL;
        headangleOL(1,cellNum) = ang_hdOL;
        headangle1OL(1,cellNum) = ang_hd1OL;
        headangle2OL(1,cellNum) = ang_hd2OL;
        
        ratemapsOR(:,:,cellNum) = rmOR;        
        ratemaps1OR(:,:,cellNum) = rm1OR;        
        ratemaps2OR(:,:,cellNum) = rm2OR;
        headangleOR(1,cellNum) = ang_hdOR;
        headangle1OR(1,cellNum) = ang_hd1OR;
        headangle2OR(1,cellNum) = ang_hd2OR;
        
        ratemapsW(:,:,cellNum) = rmW;        
        ratemaps1W(:,:,cellNum) = rm1W;        
        ratemaps2W(:,:,cellNum) = rm2W;
        headangleW(1,cellNum) = ang_hdW;
        headangle1W(1,cellNum) = ang_hd1W;
        headangle2W(1,cellNum) = ang_hd2W;
        
    end
end
%output Variable
out.time = timevar;
out.mrallO = mrlittO;
out.mrallO1 = mrlittO1;
out.mrallO2 = mrlittO2;
out.mrallOL = mrlittOL;
out.mrallO1L = mrlittO1L;
out.mrallO2L = mrlittO2L;
out.mrallOR = mrlittOR;
out.mrallO1R = mrlittO1R;
out.mrallO2R = mrlittO2R;
out.percentil99thO = prctile(out.mrallO,99);
out.correlationO = corrittO;
out.correlationOL = corrittOL;
out.correlationOR = corrittOR;
out.ratemapO = ratemapsO;
out.ratemap1O = ratemaps1O;
out.ratemap2O = ratemaps2O;
out.headangleO = headangleO;
out.headangle1O = headangle1O;
out.headangle2O = headangle2O;
out.ratemapOL = ratemapsOL;
out.ratemap1OL = ratemaps1OL;
out.ratemap2OL = ratemaps2OL;
out.headangleOL = headangleOL;
out.headangle1OL = headangle1OL;
out.headangle2OL = headangle2OL;
out.ratemapOR = ratemapsOR;
out.ratemap1OR = ratemaps1OR;
out.ratemap2OR = ratemaps2OR;
out.headangleOR = headangleOR;
out.headangle1OR = headangle1OR;
out.headangle2OR = headangle2OR;
out.mrallW = mrlittW;
out.mrallW1 = mrlittW1;
out.mrallW2 = mrlittW2;
out.percentil99thW = prctile(out.mrallW,99);
out.correlationW = corrittW;
out.ratemapW = ratemapsW;
out.ratemap1W = ratemaps1W;
out.ratemap2W = ratemaps2W;
out.headangleW = headangleW;
out.headangle1W = headangle1W;
out.headangle2W = headangle2W;
out.dimX = dimX;
out.dimY = dimY;
out.QPOL = QPOL;
out.QPOR = QPOR;
out.QPW = QPW;
out.dir = pwd;
save([pwd,'/',name,'/','EBCstats.mat'],'out');

save([pwd,'/',name,'/','disMaps.mat'],'disOL','disOR','disW','disO');
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