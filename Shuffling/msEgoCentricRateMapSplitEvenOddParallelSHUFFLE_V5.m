function [out, detailed]= msEgoCentricRateMapSplitEvenOddParallelSHUFFLE_V5(ms,HD,tracking, frameMap, dimX, dimY, deconvolve, QP, DistBinSize,DegBinSize)
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

name = ['EBCevenOddSplitParallelDeconvolved_CircSHUFFLEDonut3_100ittD',num2str(round(DistBinSize)),'A',num2str(DegBinSize),'.mat'];

if dimX > dimY
    FOVsize = round(dimY/2);
else
    FOVsize = round(dimX/2);
end

if deconvolve
    ms.FiltTraces = ms.deconvolvedSig;
end

%% Get behavior information
degBins = (-180:DegBinSize:179);                                  %Angle bins for EBC metric polar plot
degBins = degBins';                                     %reorient bins
degBins = deg2rad(degBins);                             %Convert to radians
distanceBins = 0:DistBinSize:FOVsize;                   %set to look at half the length of the field which in our case is ~38cm (37.5 rounded up)
fps = 30;                                               %Frames per second
spf = 1/fps;                                            %Seconds per frame
ms.timestamp = frameMap.*spf;                           %time stamp in seconds
minDist = [1:2];
cutoff = 300;
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

if isempty(QP)
    QP = findEdges(tracking);       %finds the edge of the graph and mark the boundaries
end
pixX = dimX/(max(QP(:,1)) - min(QP(:,1)));
pixY = dimY/(max(QP(:,2)) - min(QP(:,2)));

%% Calculate distances
degSamp = 1;                                                            %angle resolution
[dis, ex, ey] = subfunc(tracking(frameMap(:,1),1),tracking(frameMap(:,1),2),HD(frameMap), QP, degSamp);   %calls funtion to bring back wall distances when neuron fired
dis_raw = dis;
dis = fillmissing(dis,'pchip',2);                                          %interpolates missing values
[a1 a2] = find(dis<0);
% dis(a1,:) = nan;
dis = dis*pixX;                                                             %Converts boundary distances from pixels to cm.
dis = circshift(dis,90,2);                                                  %shifts values by 90 degrees

if frameMap(length(frameMap)) < length(ms.FiltTraces)
    ms.FiltTraces = ms.FiltTraces(frameMap);   %Sink trace if not already sinked
end

firing = ms.FiltTraces;
parfor itteration = 1 : 100
    shuffledFiring = CShuffle(firing);                                  %Shuffle firing peaks 
    distanceBinspar = distanceBins;
    mrlitt = zeros(length(firing(1,:)),1);
    mrlitt1 = zeros(length(firing(1,:)),1);
    mrlitt2 = zeros(length(firing(1,:)),1);
    corritt = zeros(length(firing(1,:)),1);
    %Loop through every cell, extract and analize firing instances and boundary locations
    for cellNum = 1 : length(firing(1,:))           
        fire = shuffledFiring(:,cellNum);          
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
                    dis2 = dis((1800*(i-1)+1:1800*i),:);
                elseif i < mins
                    fire2(end+1:end + 1800) = fire(1800*(i-1)+1:1800*i,1);
                    iTime2(end+1:end + 1800) =(1800*(i-1)+1:1800*i);
                    dis2(end+1:end + 1800,:) = dis(1800*(i-1)+1:1800*i,:);
                else
                    fire2(end+1:end+(length(fire)-((mins-1)*1800))) = fire((1800*(i-1)+1:end),1);
                    iTime2(end+1:end+(length(fire)-((mins-1)*1800))) = (1800*(i-1)+1:length(fire));
                    dis2(end+1:end+(length(fire)-((mins-1)*1800)),:) = dis((1800*(i-1)+1:end),:);
                end
            else
                if i == 1
                    fire1 = fire(1:1800*i,1);
                    iTime1 = (1:1800*i);
                    dis1 = dis(1:1800,:);
                elseif i < mins
                    fire1(end+1:end+1800) = fire(((1800*i)-1799:1800*i),1);
                    iTime1(end+1:end+1800)= ((1800*i)-1799:1800*i);
                    dis1(end+1:end+1800,:) = dis((1800*i)-1799:1800*i,:);
                else
                    fire1(end+1:end+(length(fire)-(mins-1)*1800)) = fire(((1800*i)-1799:end),1);
                    iTime1(end+1:end+(length(fire)-(mins-1)*1800)) = ((1800*i)-1799:length(fire));
                    dis1(end+1:end+(length(fire)-(mins-1)*1800),:) = dis((1800*i)-1799:length(fire),:);
                end
            end
        end
        
        ifire = find(fire);                                     %Find indices for all non-zero values
        ifire1 = find(fire1);                                     %Find indices for all non-zero values
        ifire2 = find(fire2);                                     %Find indices for all non-zero values
        
        if(~isempty(ifire) && ~isempty(ifire1) && ~isempty(ifire2))
            for j = 1 : length(ifire)
                %full run
                ind_fire(j) = ifire(j);                      %Add firing index to ms struct                
                %first half
                if j < length(ifire1)
                    ind_fire1(j) = ifire1(j);                      %Add firing index to ms struct                    
                end
                %second half
                if j < length(ifire2)
                    ind_fire2(j) = ifire2(j);                      %Add firing index to ms struct                    
                end
            end
            
            %% Calculate raw maps:
            thetaBins = deg2rad(linspace(-180,180,size(dis,2)));                    %angle bins
            thetaBins3d = deg2rad(linspace(-180,180,round(size(dis,2)/DegBinSize)));                    %angle bins
            occ = NaN(length(thetaBins), length(distanceBinspar));                     %wall occupancy bins
            nspk = occ;                                                             %Number of spikes bins
            nspk1 = nspk;
            nspk2 = nspk;
            occ1 = occ;
            occ2 = occ;
            distanceBinspar(end+1) = Inf;                                              %Adds an Infinity value at the end of the bins as safety procaution/break point
            ci = ind_fire(:);                                            %firing instances of the cell
            ci1 = ind_fire1(:);                                            %firing instances of the cell
            ci2 = ind_fire2(:);                                            %firing instances of the cell
            for i = 1:length(thetaBins)
                t = dis(:,i); %boundary distance for a particular bin
                t1 = dis1(:,i);
                t2 = dis2(:,i);
                for k = 1:length(distanceBinspar)-1
                    inds = t>=distanceBinspar(k) & t<distanceBinspar(k+1);                %filter through the boundary distances
                    occ(i,k) = sum(inds);                                          %Wall occupancy definition
                    inds = find(inds);                                              %find all non-zero boundary distances indices
                    nspk(i,k) = sum(fire(intersect(inds,ci)));                      %Number of spike instances definition
                    %first half
                    inds1 = t1>=distanceBinspar(k) & t1<distanceBinspar(k+1);                %filter through the boundary distances
                    occ1(i,k) = sum(inds1);                                          %Wall occupancy definition
                    inds1 = find(inds1);                                              %find all non-zero boundary distances indices
                    nspk1(i,k) = sum(fire1(intersect(inds1,ci1)));                      %Number of spike instances definition
                    %second half
                    inds2 = t2>=distanceBinspar(k) & t2<distanceBinspar(k+1);                %filter through the boundary distances
                    occ2(i,k) = sum(inds2);                                          %Wall occupancy definition
                    inds2 = find(inds2);                                              %find all non-zero boundary distances indices
                    nspk2(i,k) = sum(fire2(intersect(inds2,ci2)));                      %Number of spike instances definition
                end
            end
            occ3d = zeros(360/DegBinSize,length(occ(1,:)));
            occ13d = zeros(360/DegBinSize,length(occ(1,:)));
            occ23d = zeros(360/DegBinSize,length(occ(1,:)));
            
            nspk3d = zeros(360/DegBinSize,length(occ(1,:)));
            nspk13d = zeros(360/DegBinSize,length(occ(1,:)));
            nspk23d = zeros(360/DegBinSize,length(occ(1,:)));
            
            for i = 1 : round(length(thetaBins)/DegBinSize)
                octemp = sum(occ(i*DegBinSize-(DegBinSize-1):i*DegBinSize,:),1);
                octemp1 = sum(occ1(i*DegBinSize-(DegBinSize-1):i*DegBinSize,:),1);
                octemp2 = sum(occ2(i*DegBinSize-(DegBinSize-1):i*DegBinSize,:),1);
                
                occ3d(i,:) = octemp;
                occ13d(i,:) = octemp1;
                occ23d(i,:) = octemp2;
                
                nspktemp = sum(nspk(i*DegBinSize-(DegBinSize-1):i*DegBinSize,:),1);
                nspktemp1= sum(nspk1(i*DegBinSize-(DegBinSize-1):i*DegBinSize,:),1);
                nspktemp2 = sum(nspk2(i*DegBinSize-(DegBinSize-1):i*DegBinSize,:),1);
                
                nspk3d(i,:) = nspktemp;
                nspk13d(i,:) = nspktemp1;
                nspk23d(i,:) = nspktemp2;
            end
            
            occ = occ3d;
            occ1 = occ13d;
            occ2 = occ23d;
            
            nspk = nspk3d;
            nspk1 = nspk13d;
            nspk2 = nspk23d;
            
            distanceBinspar = distanceBinspar(1:end-1);                                   %itteration through bins
            if any(nspk(:)>0)
                
                % bring back to original dims
                occ = occ(:,1:end-1); occ=occ';
                cutout = find(occ<cutoff);
                occ(cutout) = 0;
                nspk = nspk(:,1:end-1); nspk=nspk';
                %first half
                occ1 = occ1(:,1:end-1); occ1=occ1';                
                occ1(cutout) = 0;
                nspk1 = nspk1(:,1:end-1); nspk1=nspk1';
                %second half
                occ2 = occ2(:,1:end-1); occ2=occ2';                
                occ2(cutout) = 0;
                nspk2 = nspk2(:,1:end-1); nspk2=nspk2';
                
                rm = (nspk./occ) * fps;
                rm(find(isnan(rm))) = min(rm(:));
                rm(find(isinf(rm))) = min(rm(:));
                rm = rm - min(rm(:));
                
                rm1 = (nspk1./occ1) * fps;
                rm1(find(isnan(rm1))) = min(rm1(:));
                rm1(find(isinf(rm1))) = min(rm1(:));
                rm1 = rm1 - min(rm1(:));
                
                rm2 = (nspk2./occ2) * fps;
                rm2(find(isnan(rm2))) = min(rm2(:));
                rm2(find(isinf(rm2))) = min(rm2(:));
                rm2 = rm2 - min(rm2(:));
                
                %% Smoothing
                %ratemap
                %full run                
                nd = numel(thetaBins3d);
                rm = [rm rm rm];
                rm = CMBHOME.Utils.SmoothMat(rm,smooth(1:2),smooth(3));   % Smooth it
                rm = rm(:,nd+1:2*nd); % bring it back
                %first half
                rm1 = [rm1 rm1 rm1];
                rm1 = CMBHOME.Utils.SmoothMat(rm1,smooth(1:2),smooth(3));   % Smooth it
                rm1 = rm1(:,nd+1:2*nd); % bring it back
                %second half
                rm2 = [rm2 rm2 rm2];
                rm2 = CMBHOME.Utils.SmoothMat(rm2,smooth(1:2),smooth(3));   % Smooth it
                rm2 = rm2(:,nd+1:2*nd); % bring it back
                
                rm = fliplr(rm);
                rm(minDist,:) = 0;
                rm1 = fliplr(rm1);
                rm1(minDist,:) = 0;
                rm2 = fliplr(rm2);
                rm2(minDist,:) = 0;
                corrpar = corr2(rm1,rm2);
                
                %% EBC METRIC
                %Full run
                metric = mean(rm,1)';
                if ~deconvolve
                    metric = metric - min(metric);
                end
                
                xs = metric(1:end-1).*cos(degBins(1:end-1)); % average
                ys = metric(1:end-1).*sin(degBins(1:end-1));
                
                ang_hd = atan2(mean(ys),mean(xs)); % mean direction
                
                mr = (cos(ang_hd)*sum(xs) + sin(ang_hd)*sum(ys)) / sum(metric(1:end-1)); % mean resultant length
                
                %FIRST HALF
                %EBC Metric
                metric1 = mean(rm1,1)';
                if ~deconvolve
                    metric1 = metric1 - min(metric1);
                end
                
                xs = metric1(1:end-1).*cos(degBins(1:end-1)); % average
                ys = metric1(1:end-1).*sin(degBins(1:end-1));
                
                ang_hd1 = atan2(mean(ys),mean(xs)); % mean direction
                
                mr1 = (cos(ang_hd1)*sum(xs) + sin(ang_hd1)*sum(ys)) / sum(metric1(1:end-1)); % mean resultant length
                
                %SECOND HALF
                %             %EBC Metric
                metric2 = mean(rm2,1)';
                if ~deconvolve
                    metric2 = metric2 - min(metric2);
                end
                
                xs = metric2(1:end-1).*cos(degBins(1:end-1)); % average
                ys = metric2(1:end-1).*sin(degBins(1:end-1));
                
                ang_hd2 = atan2(mean(ys),mean(xs)); % mean direction
                
                mr2 = (cos(ang_hd2)*sum(xs) + sin(ang_hd2)*sum(ys)) / sum(metric2(1:end-1)); % mean resultant length
                
                %Save Results
                ind_fire = []; %Indices of neuron activity/firing                      
                ind_fire1 = [];                
                ind_fire2 =[];      
                
                mrlitt(cellNum) = mr; %save the MRL
                mrlitt1(cellNum) =mr1;
                mrlitt2(cellNum) =mr2;
                corritt(cellNum) = corrpar;
            end
        end        
    end
    MRLsave(:,itteration) = mrlitt;
    MRLsave1(:,itteration) = mrlitt1;
    MRLsave2(:,itteration) = mrlitt2;
    CORRsave(:,itteration) = corritt;
end

if find(MRLsave>0.6)
    beep
    pause
end
out.mrall = MRLsave(:);
out.mrall1 = MRLsave1(:);
out.mrall2 = MRLsave2(:);
out.percentil99th = prctile(out.mrall,99);
out.firing = firing;
out.correlationEO = CORRsave(:);
out.dimX = dimX;
out.dimY = dimY;
out.QP = QP;

save(name,'out');
figure
histogram(out.mrall)
title('MRL distribution at 0.1 treshold')
ylabel('Number of Cells')
xlabel('Mean Resultant Length')
savefig('Shuffled_MRL_Hist_Circ_CorrDonut100.fig')
end

%% Subfunctions

%This function calculates the distance from the animal to boundaries of the environment at each behavioral data point.
%The distance calculation has to be done for all orientations around the animal centered on the animal’s
%current heading direction. That is to say that the animal’s current heading is always 0° and the distance
%to the boundaries is calculated for each of the 360 one-degree bins around the animal.
function [dis, ex, ey] = subfunc(rx,ry,hd, QP, degSamp)

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