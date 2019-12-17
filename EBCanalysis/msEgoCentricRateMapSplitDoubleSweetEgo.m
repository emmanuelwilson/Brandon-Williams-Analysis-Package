function [out, detailed]= msEgoCentricRateMapSplitDoubleSweetEgo(ms,HD,tracking, frameMap, pixX, pixY, binarize, varargin)
%%Egocentric Boundary Cell Rate Map function,boundary location polar plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function will take your data and analize it in a way to facilitate  %
%Egocentric Boundary Cell(EBC) identification. Inputs as marked above     %
%are:["ms.mat" file containing all of the physiology data, the Head       %
%Direction matrix, Head position matrix, frameMap matrix, user defined    %
%threshold (typically 0.1), x-axis pixel to cm conversion factor, and     %
%y-axis pixel to cm factor(may vary depending on video quality). varargin %
%can be ignored.                                                          %
%This function will create a new folder within you directory called "EBC  %
%results" and save all analysis figures as numbered JPG pictures in said  %
%folder.                                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author Emmanuel Wilson, modified from Jake Hinman
warning('off','stats:glmfit:IterationLimit');
warning('off','curvefit:fit:noStartPoint');
warning('off','curvefit:fittype:sethandles:WeibullxMustBePositive');
warning('off','stats:glmfit:BadScaling');
warning('off','MATLAB:nargchk:deprecated');

FOVsize = 30;
name = 'EBCresultsSweetMilk_Wall&Object_split';
%% Setup & Parse
p = inputParser;
p.addParameter('videoSamp', 1);                    % calculate every X frames of video
p.addParameter('degSamp', 1);                      % Degree bins
p.addParameter('heading', 1);                       % use heading (0--> Head direction)
p.addParameter('ifVideo', 0);                       % Play the video?
p.addParameter('labels', 1);                        % Plot forward?
p.addParameter('distanceBins', 0:1:FOVsize);           % How far to look (cm)
p.addParameter('boundaryMode', 1);                  % 0-> autolines, 1-> click, mat->useit
p.addParameter('ifLine', 0);
p.addParameter('figures',[0 0 0 1 1 1 1 1]);
p.addParameter('mergeFigures',1);                   % if 1 then puts in single figure
p.addParameter('smooth', [5 5 5])
p.parse(varargin{:});

%% Get behavior information
ratemapsO = zeros(FOVsize,361,length(ms.FiltTraces(1,:)));        %Probability ratemap values
ratemapsO1 = zeros(FOVsize,361,length(ms.FiltTraces(1,:)));        %Probability ratemap values
ratemapsO2 = zeros(FOVsize,361,length(ms.FiltTraces(1,:)));        %Probability ratemap values

ratemaps = zeros(FOVsize,361,length(ms.FiltTraces(1,:)));        %Probability ratemap values
ratemapsW1 = zeros(FOVsize,361,length(ms.FiltTraces(1,:)));        %Probability ratemap values
ratemapsW2 = zeros(FOVsize,361,length(ms.FiltTraces(1,:)));        %Probability ratemap values

ms.ind_fire = NaN(ms.numFrames,length(ms.FiltTraces(1,:))); %Indices of neuron activity/firing
ind_fire1 = NaN(ms.numFrames,length(ms.FiltTraces(1,:)));
ind_fire2 = NaN(ms.numFrames,length(ms.FiltTraces(1,:)));
ms.cell_x = NaN(ms.numFrames,length(ms.FiltTraces(1,:)));   %X coordinate of locations where cell fired
ms.cell_y = NaN(ms.numFrames,length(ms.FiltTraces(1,:)));   %Y cooridnates of locations where cell fired
ms.cell_timeo = NaN(ms.numFrames,length(ms.FiltTraces(1,:)));%Time at when cell fired
degBins = (-180:2:179);                                  %Angle bins for EBC metric polar plot
degBins = degBins';                                     %reorient bins
degBins = deg2rad(degBins);                             %Convert to radians
degBins1 = (-180:2:179);                                  %Angle bins for EBC metric polar plot
degBins1 = degBins;                                     %reorient bins
% degBins1 = deg2rad(degBins);                             %Convert to radians
degBins2 = (-180:2:179);                                  %Angle bins for EBC metric polar plot
degBins2 = degBins;                                     %reorient bins
% degBins2 = deg2rad(degBins);                             %Convert to radians
freqFire = zeros(1,length(ms.FiltTraces(1,:)));             %Firing frequency for each cell
mrall = freqFire;                                       %MRL for each cell
freqMax = freqFire;                                     %Max frequency of each cell
ms.HDfiring = ms.ind_fire;                              %Indices of neuron activity for head direction
distanceBins = 0:1:FOVsize;                                  %set to look at half the length of the field which in our case is ~38cm (37.5 rounded up)
mkdir(name)             %Create new folder within current directory
mrtotO = 0;                                              %total MRL
mrtotO1 = 0;
mrtotO2 = 0;
mrtotW = 0;                                              %total MRL
mrtotW1 = 0;
mrtotW2 = 0;
counter = 0 ;                                           %counter
counter1 = 0;
counter2 = 0;
corrO = 0;
fps = 30;                                               %Frames per second
spf = 1/fps;                                            %Seconds per frame
ms.timestamp = frameMap.*spf;                           %time stamp in seconds
county = 0;

if length(frameMap)> length(tracking(:,1))
    i = 0;
    test = frameMap(1:find(frameMap == length(tracking(:,1))));
    while isempty(test)
        test = frameMap(1:find(frameMap == length(tracking(:,1))-i));
        i = i +1;
    end
    frameMap = test;
    fprintf('FrameMap is larger than the behav')
    
end
%% Get structure of environment
%Identify where the bounds of the environment are located. Through a subfunction that allows the user to click
%on a plot of the positional data to indicate where the corners of the environment are located.

QPO = findEdges(tracking);       %finds the edge of the graph and mark the boundaries
%% Calculate distances
degSamp = 1;                                                            %angle resolution
[disO, ex, ey] = subfuncO(tracking(frameMap(:,1),1),tracking(frameMap(:,1),2),HD(frameMap), QPO, degSamp);   %calls funtion to bring back wall distances when neuron fired
dis_rawO = disO;
disO = disO*pixX;                                                             %Converts boundary distances from pixels to cm.
disO = circshift(disO,90,2);                                                  %shifts values by 90 degrees

QPW = findEdges(tracking);       %finds the edge of the graph and mark the boundaries
%% Calculate distances
degSamp = 1;                                                            %angle resolution
[disW, ex, ey] = subfuncW(tracking(frameMap(:,1),1),tracking(frameMap(:,1),2),HD(frameMap), QPW, degSamp);   %calls funtion to bring back wall distances when neuron fired
dis_rawW = disW;
disW = fillmissing(disW,'pchip',2);                                          %interpolates missing values
disW = disW*pixX;                                                             %Converts boundary distances from pixels to cm.
disW = circshift(disW,90,2);                                                  %shifts values by 90 degrees

firePeaks = Binarize(ms);
%Loop through every cell, extract and analize firing instances and boundary locations
for cellNum = 1 : length(ms.FiltTraces(1,:))
    if binarize == 1
        firing = firePeaks.binarizedTraces(:,cellNum);%ms.FiltTraces(:,cellNum);%ms.firing(:,cellNum);                          %Extract firing trace
        firing = circshift(firing,-6,1);
        firing(end-6 : end, :) = 0;        
    else
        firing = ms.FiltTraces(:,cellNum);
    end
    if length(frameMap) < length(firing)
        firing = firing(1:length(frameMap));
    end
    fire = firing;                                  %Duplicate
    fire1 = firing(1:(round(length(fire)/2)),1);
    fire2 = firing((round(length(fire)/2))+1:end,1);
    if min(fire)<0
        fire = fire + abs(min(fire)); 
        fire1 = fire1 + abs(min(fire));
        fire2 = fire2 + abs(min(fire));
    else
        fire = fire - min(fire);
        fire1 = fire1 - min(fire);
        fire2 = fire2 - min(fire);
    end
    ifire = find(fire);                                     %Find indices for all non-zero values 
    ifire1 = find(fire1);
    ifire2 = find(fire2)+length(fire1);
    if(~isempty(ifire) && ~isempty(ifire1) && ~isempty(ifire2))
        for j = 1 : length(ifire)
            %full run
            ms.ind_fire(j,cellNum) = ifire(j);                      %Add firing index to ms struct
            ms.cell_x(j,cellNum) = tracking((frameMap(ifire(j))));  %X postion of mouse during firing at sinked time
            ms.cell_y(j,cellNum) = tracking((frameMap(ifire(j))),2);%Y position of mouse during firing at sinked time
            ms.cell_time(j,cellNum) = ms.timestamp(ifire(j));       %Physiological time of firing
            ms.HDfiring(j,cellNum) = HD(frameMap(ifire(j)));        %Head Direction of mouse at time of neural firing
            %first half
            if j <= length(ifire1)
                ind_fire1(j,cellNum) = ifire1(j);                      %Add firing index to ms struct
                cell_x1(j,1) = tracking((frameMap(ifire1(j))));  %X postion of mouse during firing at sinked time
                cell_y1(j,1) = tracking((frameMap(ifire1(j))),2);%Y position of mouse during firing at sinked time
                cell_time1(j,cellNum) = ms.timestamp(ifire1(j));       %Physiological time of firing
                HDfiring1(j,cellNum) = HD(frameMap(ifire1(j)));        %Head Direction of mouse at time of neural firing
            end
            %second half
            if j > length(ifire1)
                tO = j-length(ifire1);
                ind_fire2(tO,cellNum) = ifire2(tO);                      %Add firing index to ms struct
                cell_x2(tO,1) = tracking((frameMap(ifire2(tO))));  %X postion of mouse during firing at sinked time
                cell_y2(tO,1) = tracking((frameMap(ifire2(tO))),2);%Y position of mouse during firing at sinked time
                cell_time2(tO,cellNum) = ms.timestamp(ifire2(tO));       %Physiological time of firing
                HDfiring2(tO,cellNum) = HD(frameMap(ifire2(tO)));        %Head Direction of mouse at time of neural firing
            end
        end     
        disO1 = disO(1:length(fire1),:);
        disO2 = disO(length(fire1)+1:end,:);
        disW1 = disW(1:length(fire1),:);
        disW2 = disW(length(fire1)+1:end,:);
        
        %% Calculate raw maps:
        thetaBins = deg2rad(linspace(-180,180,size(disO,2)));                    %angle bins
        occO = NaN(length(thetaBins), length(distanceBins));                     %wall occupancy bins
        nspkO = occO;                                                             %Number of spikes bins
        nspkO1 = nspkO;
        nspkO2 = nspkO;
        occO1 = occO;
        occO2 = occO;
        
        occW = NaN(length(thetaBins), length(distanceBins));                     %wall occupancy bins
        nspkW = occW;                                                             %Number of spikes bins
        nspkW1 = nspkW;
        nspkW2 = nspkW;
        occW1 = occW;
        occW2 = occW;
        
        thetaBins1 = thetaBins;
        thetaBins2 = thetaBins;
        distanceBins(end+1) = Inf;                                              %Adds an Infinity value at the end of the bins as safety procaution/break point
        distanceBins1 = distanceBins;
        distanceBins2 = distanceBins;
        ci = ms.ind_fire(:,cellNum);                                            %firing instances of the cell
        ci1 = ind_fire1(:,cellNum);                                            %firing instances of the cell
        ci2 = ind_fire2(:,cellNum);                                            %firing instances of the cell        
        for i = 1:length(thetaBins)
            tO = disO(:,i); %boundary distance for a particular bin
            tO1 = disO1(:,i);
            tO2 = disO2(:,i);
            
            t = disW(:,i); %boundary distance for a particular bin
            t1 = disW1(:,i);
            t2 = disW2(:,i);
            for k = 1:length(distanceBins)-1
                %Object
                inds = tO>=distanceBins(k) & tO<distanceBins(k+1);                %filter through the boundary distances
                occO(i,k) = sum(inds);                                          %Wall occupancy definition
                inds = find(inds);                                              %find all non-zero boundary distances indices
                nspkO(i,k) = sum(fire(intersect(inds,ci)));                      %Number of spike instances definition
                %first half
                inds1 = tO1>=distanceBins1(k) & tO1<distanceBins1(k+1);                %filter through the boundary distances
                occO1(i,k) = sum(inds1);                                          %Wall occupancy definition
                inds1 = find(inds1);                                              %find all non-zero boundary distances indices
                nspkO1(i,k) = sum(fire1(intersect(inds1,ci1)));                      %Number of spike instances definition
                %second half
                inds2 = tO2>=distanceBins2(k) & tO2<distanceBins2(k+1);                %filter through the boundary distances
                occO2(i,k) = sum(inds2);                                          %Wall occupancy definition
                inds2 = find(inds2);                                              %find all non-zero boundary distances indices
                inds2 = inds2 + length(fire1);
                nspkO2(i,k) = sum(fire2(intersect(inds2,ci2)-length(fire1)));                      %Number of spike instances definition
                
                %Wall
                inds = t>=distanceBins(k) & t<distanceBins(k+1);                %filter through the boundary distances
                occW(i,k) = sum(inds);                                          %Wall occupancy definition
                inds = find(inds);                                              %find all non-zero boundary distances indices
                nspkW(i,k) = sum(fire(intersect(inds,ci)));                      %Number of spike instances definition
                %first half
                inds1 = t1>=distanceBins1(k) & t1<distanceBins1(k+1);                %filter through the boundary distances
                occW1(i,k) = sum(inds1);                                          %Wall occupancy definition
                inds1 = find(inds1);                                              %find all non-zero boundary distances indices
                nspkW1(i,k) = sum(fire1(intersect(inds1,ci1)));                      %Number of spike instances definition
                %second half
                inds2 = t2>=distanceBins2(k) & t2<distanceBins2(k+1);                %filter through the boundary distances
                occW2(i,k) = sum(inds2);                                          %Wall occupancy definition
                inds2 = find(inds2);                                              %find all non-zero boundary distances indices
                inds2 = inds2 + length(fire1);
                nspkW2(i,k) = sum(fire2(intersect(inds2,ci2)-length(fire1)));                      %Number of spike instances definition
            end
        end
        distanceBins = distanceBins(1:end-1);                                   %itteration through bins
        distanceBins1 = distanceBins1(1:end-1);
        distanceBins2 = distanceBins2(1:end-1);
        if any(nspkO(:)>0) && any(nspkW(:)>0)
            counter = counter + 1;
            counter1 = counter1 +1;
            counter2 = counter2 +1;
            
            %Object
            % bring back to original dims
            occO = occO(:,1:end-1); occO=occO';
            nspkO = nspkO(:,1:end-1); nspkO=nspkO';
            %first half
            occO1 = occO1(:,1:end-1); occO1=occO1';
            nspkO1 = nspkO1(:,1:end-1); nspkO1=nspkO1';
            %second half
            occO2 = occO2(:,1:end-1); occO2=occO2';
            nspkO2 = nspkO2(:,1:end-1); nspkO2=nspkO2';
            
            rmO = (nspkO./occO) * fps;            
            rmO(find(isnan(rmO))) = min(rmO(:));
            rmO(find(isinf(rmO))) = min(rmO(:));
            rmO = rmO - min(rmO(:));
                        
            rmO1 = (nspkO1./occO1) * fps;            
            rmO1(find(isnan(rmO1))) = min(rmO1(:));
            rmO1(find(isinf(rmO1))) = min(rmO1(:));
            rmO1 = rmO1 - min(rmO1(:));
            
            rmO2 = (nspkO2./occO2) * fps;            
            rmO2(find(isnan(rmO2))) = min(rmO2(:));
            rmO2(find(isinf(rmO2))) = min(rmO2(:));
            rmO2 = rmO2 - min(rmO2(:));
            
            %Wall
            % bring back to original dims
            occW = occW(:,1:end-1); occW=occW';
            cutout = find(occW<25);
            occW(cutout) = 0;
            nspkW = nspkW(:,1:end-1); nspkW=nspkW';
            %first half
            occW1 = occW1(:,1:end-1); occW1=occW1';
            cutout1 = find(occW<25);
            occW1(cutout1) = 0;
            nspkW1 = nspkW1(:,1:end-1); nspkW1=nspkW1';
            %second half
            occW2 = occW2(:,1:end-1); occW2=occW2';
            cutout2 = find(occW<25);
            occW2(cutout2) = 0;
            nspkW2 = nspkW2(:,1:end-1); nspkW2=nspkW2';
            
            rmW = (nspkW./occW) * fps;            
            rmW(find(isnan(rmW))) = min(rmW(:));
            rmW(find(isinf(rmW))) = min(rmW(:));
            rmW = rmW - min(rmW(:));
                        
            rmW1 = (nspkW1./occW1) * fps;            
            rmW1(find(isnan(rmW1))) = min(rmW1(:));
            rmW1(find(isinf(rmW1))) = min(rmW1(:));
            rmW1 = rmW1 - min(rmW1(:));
            
            rmW2 = (nspkW2./occW2) * fps;            
            rmW2(find(isnan(rmW2))) = min(rmW2(:));
            rmW2(find(isinf(rmW2))) = min(rmW2(:));
            rmW2 = rmW2 - min(rmW2(:));
            
            %% Smoothing
            %Object
            %occupancy
            %full run
            occO = [occO occO occO];
            nd = numel(thetaBins);
            occO = CMBHOME.Utils.SmoothMat(occO, smooth(1:2), smooth(3));
            occO = occO(:, nd+1:2*nd);
            %first half
            occO1 = [occO1 occO1 occO1];
            nd1 = numel(thetaBins1);
            occO1 = CMBHOME.Utils.SmoothMat(occO1, smooth(1:2), smooth(3));
            occO1 = occO1(:, nd1+1:2*nd1);
            %second half
            occO2 = [occO2 occO2 occO2];
            nd2 = numel(thetaBins2);
            occO2 = CMBHOME.Utils.SmoothMat(occO2, smooth(1:2), smooth(3));
            occO2 = occO2(:, nd2+1:2*nd2);
            
            %number of spikes
            %full run
            nspkO = [nspkO nspkO nspkO];
            nspkO = CMBHOME.Utils.SmoothMat(nspkO,smooth(1:2),smooth(3));   % Smooth it
            nspkO = nspkO(:,nd+1:2*nd); % bring it back
            %first half
            nspkO1 = [nspkO1 nspkO1 nspkO1];
            nspkO1 = CMBHOME.Utils.SmoothMat(nspkO1,smooth(1:2),smooth(3));   % Smooth it
            nspkO1 = nspkO1(:,nd1+1:2*nd1); % bring it back
            %second half
            nspkO2 = [nspkO2 nspkO2 nspkO2];
            nspkO2 = CMBHOME.Utils.SmoothMat(nspkO2,smooth(1:2),smooth(3));   % Smooth it
            nspkO2 = nspkO2(:,nd2+1:2*nd2); % bring it back
            
            %ratemap
            %full run
%             rm = (nspk./occ) * fps;
            rmO = [rmO rmO rmO];
            rmO = CMBHOME.Utils.SmoothMat(rmO,smooth(1:2),smooth(3));   % Smooth it
            rmO = rmO(:,nd+1:2*nd); % bring it back
            %first half
            rmO1 = [rmO1 rmO1 rmO1];
            rmO1 = CMBHOME.Utils.SmoothMat(rmO1,smooth(1:2),smooth(3));   % Smooth it
            rmO1 = rmO1(:,nd1+1:2*nd1); % bring it back
            %second half
            rmO2 = [rmO2 rmO2 rmO2];
            rmO2 = CMBHOME.Utils.SmoothMat(rmO2,smooth(1:2),smooth(3));   % Smooth it
            rmO2 = rmO2(:,nd2+1:2*nd2); % bring it back
            
            occO = fliplr(occO);
            nspkO = fliplr(nspkO);
            rmO = fliplr(rmO);
            
            occO1 = fliplr(occO1);
            nspkO1 = fliplr(nspkO1);
            rmO1 = fliplr(rmO1);
            
            occO2 = fliplr(occO2);
            nspkO2 = fliplr(nspkO2);
            rmO2 = fliplr(rmO2);
            
            %Wall
            %occupancy
            %full run
            occW = [occW occW occW];
            nd = numel(thetaBins);
            occW = CMBHOME.Utils.SmoothMat(occW, smooth(1:2), smooth(3));
            occW = occW(:, nd+1:2*nd);
            occW(cutout) = 0;
            %first half
            occW1 = [occW1 occW1 occW1];
            nd1 = numel(thetaBins1);
            occW1 = CMBHOME.Utils.SmoothMat(occW1, smooth(1:2), smooth(3));
            occW1 = occW1(:, nd1+1:2*nd1);
            occW1(cutout) = 0;
            %second half
            occW2 = [occW2 occW2 occW2];
            nd2 = numel(thetaBins2);
            occW2 = CMBHOME.Utils.SmoothMat(occW2, smooth(1:2), smooth(3));
            occW2 = occW2(:, nd2+1:2*nd2);
            occW2(cutout) = 0;
            
            %number of spikes
            %full run
            nspkW = [nspkW nspkW nspkW];
            nspkW = CMBHOME.Utils.SmoothMat(nspkW,smooth(1:2),smooth(3));   % Smooth it
            nspkW = nspkW(:,nd+1:2*nd); % bring it back
            %first half
            nspkW1 = [nspkW1 nspkW1 nspkW1];
            nspkW1 = CMBHOME.Utils.SmoothMat(nspkW1,smooth(1:2),smooth(3));   % Smooth it
            nspkW1 = nspkW1(:,nd1+1:2*nd1); % bring it back
            %second half
            nspkW2 = [nspkW2 nspkW2 nspkW2];
            nspkW2 = CMBHOME.Utils.SmoothMat(nspkW2,smooth(1:2),smooth(3));   % Smooth it
            nspkW2 = nspkW2(:,nd2+1:2*nd2); % bring it back
            
            %ratemap
            %full run
%             rm = (nspk./occ) * fps;
            rmW = [rmW rmW rmW];
            rmW = CMBHOME.Utils.SmoothMat(rmW,smooth(1:2),smooth(3));   % Smooth it
            rmW = rmW(:,nd+1:2*nd); % bring it back
            %first half
            rmW1 = [rmW1 rmW1 rmW1];
            rmW1 = CMBHOME.Utils.SmoothMat(rmW1,smooth(1:2),smooth(3));   % Smooth it
            rmW1 = rmW1(:,nd1+1:2*nd1); % bring it back
            %second half
            rmW2 = [rmW2 rmW2 rmW2];
            rmW2 = CMBHOME.Utils.SmoothMat(rmW2,smooth(1:2),smooth(3));   % Smooth it
            rmW2 = rmW2(:,nd2+1:2*nd2); % bring it back
            
            occW = fliplr(occW);
            nspkW = fliplr(nspkW);
            rmW = fliplr(rmW);
            
            occW1 = fliplr(occW1);
            nspkW1 = fliplr(nspkW1);
            rmW1 = fliplr(rmW1);
            
            occW2 = fliplr(occW2);
            nspkW2 = fliplr(nspkW2);
            rmW2 = fliplr(rmW2);
            
            %% Plots
            %The first three plots are rectangular versions of the three color-coded plots
            %Figure 7 is part of our work in progress for estimating a cell’s preferred distance.
            %Finally, figure 8 is a trajectory plot with heading color coded spike dots
            
            figure(1);
            n=3;
            m = 5;
            c = 1;                                              
            
            % ratemap circular object
            subplot(n,m,c);c=c+1;
            % the +pi/2 brings "forwards" to "up"
            [tO2, r2] = meshgrid(wrapTo2Pi(thetaBins+pi/2), distanceBins(1:end-1));
            [x, y] = pol2cart(tO2,r2);
            h=surface(x,y, rmO); shading interp                       
            hold on
            set(gca,'XTick',[],'YTick',[])
            axis square
            axis off
            colormap(jet)
            set(gca,'YDir','Normal')
            freezeColors
            title('rm Obj')
            
             % ratemap circular wall
            subplot(n,m,c);c=c+1;
            % the +pi/2 brings "forwards" to "up"
            [t2, r2] = meshgrid(wrapTo2Pi(thetaBins+pi/2), distanceBins(1:end-1));
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
            
            %Trajectory map
            if binarize == 1
                subplot(n,m,c);c=c+1;
                edg = splitter(QPO);
                hold on
                plot(pixX*tracking(:,1),-pixY*tracking(:,2),'Color',[.7 .7 .7])
                colormap(hsv)
                xlim(pixX*[min(tracking(:,1)) max(tracking(:,1))]);ylim(pixY*[-max(tracking(:,2)) -min(tracking(:,2))])
                cx=pixX*ms.cell_x(:,cellNum);
                cy=pixY*ms.cell_y(:,cellNum);
                scatter(cx,-cy,38,ms.HDfiring(:,cellNum),'filled')
                set(gca,'YDir','Normal')
                caxis([0 360])
                title('Traj')
                axis off
                axis square
            else
                subplot(n,m,c);c=c+1;
                edg = splitter(QPO);
                hold on
                plot(pixX*tracking(:,1),-pixY*tracking(:,2),'Color',[.7 .7 .7])
                colormap(hsv)
                xlim(pixX*[min(tracking(:,1)) max(tracking(:,1))]);ylim(pixY*[-max(tracking(:,2)) -min(tracking(:,2))])
                cx=pixX*ms.cell_x(:,cellNum);
                cy=pixY*ms.cell_y(:,cellNum);
                
                for scatobj = 1 : 5
                    perc = scatobj*20;
                    scatTrace = fire;
                    scatTrace(scatTrace> prctile(scatTrace,perc) | scatTrace<= prctile(scatTrace,(perc-20)))=0;
                    scatInd = find(scatTrace);
                    s=scatter(cx(scatInd),-cy(scatInd),38,ms.HDfiring(scatInd,cellNum),'filled');
                    s.MarkerFaceAlpha = (perc-20)*0.01;
                end                   
            end
            set(gca,'YDir','Normal')    
            caxis([0 360])
            axis([min(pixX*tracking(:,1)) max(pixX*tracking(:,1)) min(-pixY*tracking(:,2)) max(-pixY*tracking(:,2))])
            title('Traj')
            axis off
            axis square
            
%             %EBC Metric
            avgcount = zeros(1,i);
            metric = zeros(1,180);
            subplot(n,m,c);c=c+1;

            r = 0;
            for it = 1 : i
                avgcount(1,it) = mean(rmO(:,it));
                if mod(it,2) == 0
                    r = it/2;
                    metric(1,r) = (avgcount(1,it-1)+avgcount(1,it))/2;
                end
            end
            metric = metric';
            if ~binarize                
                metric = metric - min(metric);
            end
            polarplot(degBins,metric)
            
            xs = metric(1:end-1).*cos(degBins(1:end-1)); % average
            ys = metric(1:end-1).*sin(degBins(1:end-1));
            
            coordlims=axis;
            
            ang_hd = atan2(mean(ys),mean(xs)); % mean direction
            
            mrO = (cos(ang_hd)*sum(xs) + sin(ang_hd)*sum(ys)) / sum(metric(1:end-1)); % mean resultant length
            
            mag_hd = sqrt(sum(ys)^2+sum(xs)^2)/sqrt(sum(abs(ys))^2+sum(abs(xs))^2)*6.28; % for visualizations sake
            
            hold on;
            polarplot([ang_hd ang_hd ],[0 mrO], 'r')
            pol = gca;
            pol.ThetaZeroLocation = 'top';
            pol.ThetaColor = 'none';
            pol.RTickLabel = [];
            hold off
            title('Wall Directionality Obj')
            stat = ['MRL: ' num2str(mrO) 'Angle : ' num2str(rad2deg(ang_hd))];
            text(0.2,coordlims(4),stat);              
            
            %EBC Metric Wall
            avgcount = zeros(1,i);
            metric = zeros(1,180);
            subplot(n,m,c);c=c+1;

            r = 0;
            for it = 1 : i
                avgcount(1,it) = mean(rmW(:,it));
                if mod(it,2) == 0
                    r = it/2;
                    metric(1,r) = (avgcount(1,it-1)+avgcount(1,it))/2;
                end
            end
            metric = metric';
            if ~binarize                
                metric = metric - min(metric);
            end
            polarplot(degBins,metric)
            
            xs = metric(1:end-1).*cos(degBins(1:end-1)); % average
            ys = metric(1:end-1).*sin(degBins(1:end-1));
            
            coordlims=axis;
            
            ang_hd = atan2(mean(ys),mean(xs)); % mean direction
            
            mrW = (cos(ang_hd)*sum(xs) + sin(ang_hd)*sum(ys)) / sum(metric(1:end-1)); % mean resultant length
            
            mag_hd = sqrt(sum(ys)^2+sum(xs)^2)/sqrt(sum(abs(ys))^2+sum(abs(xs))^2)*6.28; % for visualizations sake
            
            hold on;
            polarplot([ang_hd ang_hd ],[0 mrW], 'r')
            pol = gca;
            pol.ThetaZeroLocation = 'top';
            pol.ThetaColor = 'none';
            pol.RTickLabel = [];
            hold off
            title('Wall Directionality Wall')
            stat = ['MRL: ' num2str(mrW) 'Angle : ' num2str(rad2deg(ang_hd))];
            text(0.2,coordlims(4),stat);                          
            
            %FIRST HALF                              
            % ratemap circular Object
            subplot(n,m,c);c=c+1;
            % the +pi/2 brings "forwards" to "up"
            [tO2, r2] = meshgrid(wrapTo2Pi(thetaBins1+pi/2), distanceBins1(1:end-1));
            [x, y] = pol2cart(tO2,r2);
            h=surface(x,y, rmO1); shading interp

            hold on
            set(gca,'XTick',[],'YTick',[])
            axis square
            axis off
            colormap(jet)
            set(gca,'YDir','Normal')
            freezeColors
            title('rm Obj 1st half')
         
            % ratemap circular Wall
            subplot(n,m,c);c=c+1;
            % the +pi/2 brings "forwards" to "up"
            [t2, r2] = meshgrid(wrapTo2Pi(thetaBins1+pi/2), distanceBins1(1:end-1));
            [x, y] = pol2cart(t2,r2);
            h=surface(x,y, rmW1); shading interp
            hold on
            set(gca,'XTick',[],'YTick',[])
            axis square
            axis off
            colormap(jet)
            set(gca,'YDir','Normal')
            freezeColors
            title('rm Wall 1st half')
            
            %Trajectory map
            if binarize == 1
                subplot(n,m,c);c=c+1;
                edg = splitter(QPO);
                hold on
                plot(pixX*tracking(1:round(length(tracking)/2),1),-pixY*tracking(1:round(length(tracking)/2),2),'Color',[.7 .7 .7])
                colormap(hsv)
                xlim(pixX*[min(tracking(:,1)) max(tracking(:,1))]);ylim(pixY*[-max(tracking(:,2)) -min(tracking(:,2))])
                cx=pixX*cell_x1(:,1);
                cy=pixY*cell_y1(:,1);
                scatter(cx,-cy,38,HDfiring1(:,cellNum),'filled')
                set(gca,'YDir','Normal')
                caxis([0 360])
                title('Traj')
                axis off
                axis square
            else
                subplot(n,m,c);c=c+1;
                edg = splitter(QPO);
                hold on
                plot(pixX*tracking(:,1),-pixY*tracking(:,2),'Color',[.7 .7 .7])
                colormap(hsv)
                xlim(pixX*[min(tracking(:,1)) max(tracking(:,1))]);ylim(pixY*[-max(tracking(:,2)) -min(tracking(:,2))])
                cx=pixX*cell_x1(:,1);
                cy=pixY*cell_y1(:,1);
                
                for scatobj = 1 : 5
                    perc = scatobj*20;
                    scatTrace = fire1;
                    scatTrace(scatTrace> prctile(scatTrace,perc) | scatTrace<= prctile(scatTrace,(perc-20)))=0;
                    scatInd = find(scatTrace);
                    s=scatter(cx(scatInd),-cy(scatInd),38,HDfiring(scatInd,cellNum),'filled');
                    s.MarkerFaceAlpha = (perc-20)*0.01;
                end                   
            end
            set(gca,'YDir','Normal')    
            caxis([0 360])
            axis([min(pixX*tracking(:,1)) max(pixX*tracking(:,1)) min(-pixY*tracking(:,2)) max(-pixY*tracking(:,2))])
            title('Traj')
            axis off
            axis square
            
%             %EBC Metric Object
            avgcount1 = zeros(1,i);
            metric1 = zeros(1,180);
            subplot(n,m,c);c=c+1;

            r = 0;
            for it = 1 : i
                avgcount1(1,it) = mean(rmO1(:,it));
                if mod(it,2) == 0
                    r = it/2;
                    metric1(1,r) = (avgcount1(1,it-1)+avgcount1(1,it))/2;
                end
            end
            
            metric1 = metric1';
            if ~binarize                
                metric1 = metric1 - min(metric1);
            end
            polarplot(degBins1,metric1)
            
            xs = metric1(1:end-1).*cos(degBins1(1:end-1)); % average
            ys = metric1(1:end-1).*sin(degBins1(1:end-1));
            
            coordlims=axis;
            
            ang_hd1 = atan2(mean(ys),mean(xs)); % mean direction
            
            mrO1 = (cos(ang_hd1)*sum(xs) + sin(ang_hd1)*sum(ys)) / sum(metric1(1:end-1)); % mean resultant length
            
            mag_hd1 = sqrt(sum(ys)^2+sum(xs)^2)/sqrt(sum(abs(ys))^2+sum(abs(xs))^2)*6.28; % for visualizations sake
            
            hold on;
            polarplot([ang_hd1 ang_hd1 ],[0 mrO1], 'r')
            pol = gca;
            pol.ThetaZeroLocation = 'top';
            pol.ThetaColor = 'none';
            pol.RTickLabel = [];
            hold off
%             set(gca,'Color','none')
            title('Wall Directionality Obj: 1st half')
            stat = ['MRL: ' num2str(mrO1) 'Angle : ' num2str(rad2deg(ang_hd))];
            text(0.2,coordlims(4),stat);
            
            %EBC Metric Wall
            avgcount1 = zeros(1,i);
            metric1 = zeros(1,180);
            subplot(n,m,c);c=c+1;

            r = 0;
            for it = 1 : i
                avgcount1(1,it) = mean(rmW1(:,it));
                if mod(it,2) == 0
                    r = it/2;
                    metric1(1,r) = (avgcount1(1,it-1)+avgcount1(1,it))/2;
                end
            end
            metric1 = metric1';
            if ~binarize                
                metric1 = metric1 - min(metric1);
            end
            polarplot(degBins1,metric1)
            
            xs = metric1(1:end-1).*cos(degBins1(1:end-1)); % average
            ys = metric1(1:end-1).*sin(degBins1(1:end-1));
            
            coordlims=axis;
            
            ang_hd1 = atan2(mean(ys),mean(xs)); % mean direction
            
            mrW1 = (cos(ang_hd1)*sum(xs) + sin(ang_hd1)*sum(ys)) / sum(metric1(1:end-1)); % mean resultant length
            
            mag_hd1 = sqrt(sum(ys)^2+sum(xs)^2)/sqrt(sum(abs(ys))^2+sum(abs(xs))^2)*6.28; % for visualizations sake
            
            hold on;
            polarplot([ang_hd1 ang_hd1 ],[0 mrW1], 'r')
            pol = gca;
            pol.ThetaZeroLocation = 'top';
            pol.ThetaColor = 'none';
            pol.RTickLabel = [];
            hold off
            title('Wall Directionality Wall: 1st half')
            stat = ['MRL: ' num2str(mrW1) 'Angle : ' num2str(rad2deg(ang_hd))];
            text(0.2,coordlims(4),stat);
            
            %SECOND HALF            
            % ratemap circular Object
            subplot(n,m,c);c=c+1;
            % the +pi/2 brings "forwards" to "up"
            [tO2, r2] = meshgrid(wrapTo2Pi(thetaBins2+pi/2), distanceBins2(1:end-1));
            [x, y] = pol2cart(tO2,r2);
            h=surface(x,y, rmO2); shading interp            
            hold on
            set(gca,'XTick',[],'YTick',[])
            axis square
            axis off
            colormap(jet)
            set(gca,'YDir','Normal')
            freezeColors
            title('rm Obj 2nd half')
            
            % ratemap circular Wall
            subplot(n,m,c);c=c+1;
            % the +pi/2 brings "forwards" to "up"
            [t2, r2] = meshgrid(wrapTo2Pi(thetaBins2+pi/2), distanceBins2(1:end-1));
            [x, y] = pol2cart(t2,r2);
            h=surface(x,y, rmW2); shading interp
            hold on
            set(gca,'XTick',[],'YTick',[])
            axis square
            axis off            
            colormap(jet)
            set(gca,'YDir','Normal')
            freezeColors
            title('rm Wall 2nd half')
            
            %Trajectory map
            if binarize == 1
                subplot(n,m,c);c=c+1;
                edg = splitter(QPO);
                hold on
                plot(pixX*tracking(round(length(tracking)/2):end,1),-pixY*tracking(round(length(tracking)/2):end,2),'Color',[.7 .7 .7])
                colormap(hsv)
                xlim(pixX*[min(tracking(:,1)) max(tracking(:,1))]);ylim(pixY*[-max(tracking(:,2)) -min(tracking(:,2))])                
                cx=pixX*cell_x2(:,1);
                cy=pixY*cell_y2(:,1);
                scatter(cx,-cy,38,HDfiring2(:,cellNum),'filled')
                set(gca,'YDir','Normal')
                caxis([0 360])
                title('Traj Obj')
                axis off
                axis square
            else
                subplot(n,m,c);c=c+1;
                edg = splitter(QPO);
                hold on
                plot(pixX*tracking(round(length(tracking)/2):end,1),-pixY*tracking(round(length(tracking)/2):end,2),'Color',[.7 .7 .7])
                colormap(hsv)
                xlim(pixX*[min(tracking(:,1)) max(tracking(:,1))]);ylim(pixY*[-max(tracking(:,2)) -min(tracking(:,2))])
                cx=pixX*cell_x2(:,1);
                cy=pixY*cell_y2(:,1);                
                for scatobj = 1 : 5
                    perc = scatobj*20;
                    scatTrace = fire2;
                    scatTrace(scatTrace> prctile(scatTrace,perc) | scatTrace<= prctile(scatTrace,(perc-20)))=0;
                    scatInd = find(scatTrace);
                    s=scatter(cx(scatInd),-cy(scatInd),38,HDfiring(scatInd,cellNum),'filled');
                    s.MarkerFaceAlpha = (perc-20)*0.01;
                end                   
            end
            set(gca,'YDir','Normal')    
            caxis([0 360])
            axis([min(pixX*tracking(:,1)) max(pixX*tracking(:,1)) min(-pixY*tracking(:,2)) max(-pixY*tracking(:,2))])
            title('Traj')
            axis off
            axis square
            
%             %EBC Metric Object
            avgcount2 = zeros(1,i);
            metric2 = zeros(1,180);
            subplot(n,m,c);
            c=c+1;
            r = 0;
            for it = 1 : i
                avgcount2(1,it) = mean(rmO2(:,it));
                if mod(it,2) == 0
                    r = it/2;
                    metric2(1,r) = (avgcount2(1,it-1)+avgcount2(1,it))/2;
                end
            end
            metric2 = metric2';
            if ~binarize                
                metric2 = metric2 - min(metric2);
            end
            polarplot(degBins2,metric2)
            
            xs = metric2(1:end-1).*cos(degBins2(1:end-1)); % average
            ys = metric2(1:end-1).*sin(degBins2(1:end-1));
            
            coordlims=axis;
            
            ang_hd2 = atan2(mean(ys),mean(xs)); % mean direction
            
            mrO2 = (cos(ang_hd2)*sum(xs) + sin(ang_hd2)*sum(ys)) / sum(metric2(1:end-1)); % mean resultant length
            
            mag_hd2 = sqrt(sum(ys)^2+sum(xs)^2)/sqrt(sum(abs(ys))^2+sum(abs(xs))^2)*6.28; % for visualizations sake
            
            hold on;
            polarplot([ang_hd2 ang_hd2 ],[0 mrO2], 'r')
            pol = gca;
            pol.ThetaZeroLocation = 'top';
            pol.ThetaColor = 'none';
            pol.RTickLabel = [];
            hold off
            title('Wall Directionality Obj')
            stat = ['MRL: ' num2str(mrO2) 'Angle : ' num2str(rad2deg(ang_hd2))];
            text(0.2,coordlims(4),stat);           
            
            %EBC Metric Wall
            avgcount2 = zeros(1,i);
            metric2 = zeros(1,180);
            subplot(n,m,c);

            r = 0;
            for it = 1 : i
                avgcount2(1,it) = mean(rmW2(:,it));
                if mod(it,2) == 0
                    r = it/2;
                    metric2(1,r) = (avgcount2(1,it-1)+avgcount2(1,it))/2;
                end
            end
            metric2 = metric2';
            if ~binarize                
                metric2 = metric2 - min(metric2);
            end
            polarplot(degBins2,metric2)
            
            xs = metric2(1:end-1).*cos(degBins2(1:end-1)); % average
            ys = metric2(1:end-1).*sin(degBins2(1:end-1));
            
            coordlims=axis;
            
            ang_hd2 = atan2(mean(ys),mean(xs)); % mean direction
            
            mrW2 = (cos(ang_hd2)*sum(xs) + sin(ang_hd2)*sum(ys)) / sum(metric2(1:end-1)); % mean resultant length
            
            mag_hd2 = sqrt(sum(ys)^2+sum(xs)^2)/sqrt(sum(abs(ys))^2+sum(abs(xs))^2)*6.28; % for visualizations sake
            
            hold on;
            polarplot([ang_hd2 ang_hd2 ],[0 mrW2], 'r')
            pol = gca;
            pol.ThetaZeroLocation = 'top';
            pol.ThetaColor = 'none';
            pol.RTickLabel = [];
            hold off
            title('Wall Directionality')
            stat = ['MRL: ' num2str(mrW2) 'Angle : ' num2str(rad2deg(ang_hd2))];
            text(0.2,coordlims(4),stat);
            
            
            corrO(1,cellNum) = corr2(rmO1,rmO2);
            corrW(1,cellNum) = corr2(rmW1,rmW2);
            
            %Save Results
            saveas(gcf,[name,'/',num2str(cellNum),'EBC.jpg']); %saving figure as a picture file (.jpg) in the new folder "EBCresults"
            clf
            ms.ind_fire = NaN(ms.numFrames,length(ms.FiltTraces(1,:))); %Indices of neuron activity/firing
            ms.cell_x = NaN(ms.numFrames,length(ms.FiltTraces(1,:)));   %X coordinate of locations where cell fired
            ms.cell_y = NaN(ms.numFrames,length(ms.FiltTraces(1,:)));   %Y cooridnates of locations where cell fired
            ms.cell_time = NaN(ms.numFrames,length(ms.FiltTraces(1,:)));%Time at when cell fired
            ms.HDfiring = ms.ind_fire;                              %Indices of neuron activity for head direction 
            ind_fire1 = NaN(size(ind_fire1));
            cell_x1 = NaN;
            cell_y1 = NaN;
            cell_time1 = NaN(size(cell_time1));
            HDfiring1 = NaN;
            ind_fire2 = NaN(size(ind_fire2));
            cell_x2 = NaN;
            cell_y2 = NaN;
            cell_time2 = NaN(size(cell_time2));
            HDfiring2 = NaN;
             
            ratemapsO(:,:,cellNum) = rmO;
            mrtotO = mrtotO + mrO;
            freqFire(1,cellNum) = length(ifire);
            mrallO(1,cellNum) = mrO;
            
            ratemapsO1(:,:,cellNum) = rmO1;
            mrtotO1 = mrtotO1 + mrO1;
            freqFire1(1,cellNum) = length(ifire1);
            mrall1(1,cellNum) = mrO1; 
            
            ratemapsO2(:,:,cellNum) = rmO2;
            mrtotO2 = mrtotO2 + mrO2;
            freqFire2(1,cellNum) = length(ifire2);
            mrallO2(1,cellNum) = mrO2;
            
            ratemapsW(:,:,cellNum) = rmW;
            mrtotW = mrtotW + mrW;
            mrallW(1,cellNum) = mrW;
            
            ratemapsW1(:,:,cellNum) = rmW1;
            mrtotW1 = mrtotW1 + mrW1;
            mrallW1(1,cellNum) = mrW1; 
            
            ratemapsW2(:,:,cellNum) = rmW2;
            mrtotW2 = mrtotW2 + mrW2;
            mrallW2(1,cellNum) = mrW2;
            
        end
    end
end
out.corrObject = corrO;
out.corrWall = corrW;
out.mravgObject = mrtotO/counter;
out.mrallObject = mrallO;
out.mrallWall = mrallW;
out.freqMax = freqMax;
out.rmObject = ratemapsO;
out.rm1Object = ratemapsO1;
out.rm2Object = ratemapsO2;
out.rmWall = ratemapsW;
out.rm1Wall = ratemapsW1;
out.rm2Wall = ratemapsW2;
out.frameNum = length(frameMap(:,1));
out.QPobject = QPO;
out.QPwalls = QPW;
out.disObject = disO;
out.disWall = disW;
save('EBCstats.mat','out');
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
%         indexes = (px>bb(1,1) & px<bb(1,2) & py>bb(2,1) & py<bb(2,2)); 
        indexes = ~((px >= bb(1,1) & px <= bb(1,2)) & (py>= bb(2,1) & py <= bb(2,2)));            
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