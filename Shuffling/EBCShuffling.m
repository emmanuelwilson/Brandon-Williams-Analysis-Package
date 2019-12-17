function [out]= EBCShuffling(ms,HD,tracking, frameMap, thresh, pixX, pixY)
%%Returns Mean Resultant Length distribution of data
%INPUTS:-firing: n by m matrix where n is spike amplitude at that frame and 
%       where m represents the cell number.
%       -HD: Head direction of the mouse in degrees
%       -tracking: Mouse/Mouse head X and Y position through time
%       -frameMap: Miniscope to webcam(ie.physiology to behavriour cameras)
%       frame syncronization
%       -thresh: threshhold set by the user, any firing value below the 
%       threshold will be ignored
%       -pixX/pixY: aproximate cm/pixel converstion in the X and Y directon
%
%OUTPUT:-out:
%           +mrall: n by 1 matrix containing all shuffled MRL values
%           +percentil99th: 99th percentile value of shuffled MRL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Emmanuel Wilson

firePeaks = Binarize(ms);
firing = firePeaks.binarizedTraces;                              %Extract firing trace
firing = circshift(firing,-6,1);
firing(end-6 : end, :) = 0;

MRLsave = zeros(length(firing(1,:))*1000,1);                                 %saves all values for all MRL values
ind_fire = NaN(length(firing(:,1)),length(firing(1,:)));                    %Indices of neuron activity/firing
distanceBins = 0:1:30;                                                      %set to look at half the length of the field which in our case is ~38cm (37.5 rounded up)
degBins = (-180:2:179);                                                      %Angle bins for EBC metric polar plot
degBins = degBins';                                                         %Transposing bins
degBins = deg2rad(degBins);                                                 %converting from degrees to radians
fps = 30;                                                                   %Frames per second

%% Get structure of environment
%Identify where the bounds of the environment are located. Through a subfunction that allows the user to click
%on a plot of the positional data to indicate where the corners of the environment are located.

if length(frameMap)> length(tracking(:,1))
    frameMap = frameMap(1:find(frameMap == length(tracking(:,1))));
    fprintf('FrameMap is larger than the behav')
end
[firing, freq] = CellFiltering(firing,thresh);                              %Eliminates the bottom 5% of cells due to low firing
QP = findEdges(tracking);       %finds the edge of the graph and mark the boundaries
%% Calculate distances
degSamp = 1;                                                            %angle resolution
[dis, ex, ey] = subfunc(tracking(frameMap(:,1),1),tracking(frameMap(:,1),2),HD(frameMap), QP, degSamp);   %calls funtion to bring back wall distances when neuron fired
dis_raw = dis;
dis = fillmissing(dis,'pchip',2);                                %fills any empty spaces that were missed in subfunc
dis = dis*pixX;                                                             %Converts boundary distances from pixels to cm.
dis = circshift(dis,90,2);  

for itteration = 1 : 100                                                    %Loop the number of times you wish to shuffle
    shuffledFiring = ShufflePeaks(firing);                                  %Shuffle firing peaks 
    
    %Loop through every cell, extract and analize firing instances and boundary locations
    for cellNum = 1 : length(firing(1,:))               
        fire = shuffledFiring(:,cellNum);                                   %Look at shuffled firing
        fire(fire < thresh) = 0;                                            %apply threshold on firing
        ifire = find(fire);                                                 %Find indices for all non-zero values
        if(~isempty(ifire))
            for j = 1 : length(ifire)
                ind_fire(j,cellNum) = ifire(j);                             %Add firing index to ms struct
            end            
            
            %% Calculate raw maps:
            thetaBins = deg2rad(linspace(-180,180,size(dis,2)));            %angle bins
            occ = NaN(length(thetaBins), length(distanceBins));             %wall occupancy bins
            nspk = occ;                                                     %Number of spikes bins
            distanceBins(end+1) = Inf;                                      %Adds an Infinity value at the end of the bins as safety procaution/break point
            ci = ind_fire(:,cellNum);                                       %firing instances of the cell
            for i = 1:length(thetaBins)
                t = dis(:,i);                                               %boundary distance for a particular bin
                for k = 1:length(distanceBins)-1
                    inds = t>=distanceBins(k) & t<distanceBins(k+1);        %filter through the boundary distances
                    occ(i,k) = sum(inds);                                   %Wall occupancy definition
                    inds = find(inds);                                      %find all non-zero boundary distances indices
                    nspk(i,k) = sum(fire(intersect(inds,ci)));                 %Number of spike instances definition
                end
            end
            distanceBins = distanceBins(1:end-1);                           %Get rid of Infinity
            if any(nspk(:)>0)
                % bring back to original dims
                occ = occ(:,1:end-1); occ=occ';                             
                nspk = nspk(:,1:end-1); nspk=nspk';
                
                rm = (nspk./occ) * fps;
                rm(find(isnan(rm))) = min(rm(:));
                rm(find(isinf(rm))) = min(rm(:));
                rm = rm - min(rm(:));
                
                %% Smoothing
                occ = [occ occ occ];
                nd = numel(thetaBins);
                occ = CMBHOME.Utils.SmoothMat(occ, smooth(1:2), smooth(3));
                occ = occ(:, nd+1:2*nd);
                
                nspk = [nspk nspk nspk];
                nspk = CMBHOME.Utils.SmoothMat(nspk,smooth(1:2),smooth(3));   % Smooth it
                nspk = nspk(:,nd+1:2*nd); % bring it back
                
                
                rm = [rm rm rm];
                rm = CMBHOME.Utils.SmoothMat(rm,smooth(1:2),smooth(3));   % Smooth it
                rm = rm(:,nd+1:2*nd); % bring it back
                
                %MRL Calculations
                avgcount = zeros(1,i);                                      %average value of ratemap at each angle
                metric = zeros(1,180);                                      %average value of ratemap for every second angle
                for it = 1 : i
                    avgcount(1,it) = mean(rm(:,it));                        %average value along angle it
                    if mod(it,2) == 0                                       %for every second angle
                        r = it/2;                                           %half the angle value
                        metric(1,r) = (avgcount(1,it-1)+avgcount(1,it))/2;  %average between current angle and previous angle
                    end
                end
                metric = metric';
%                 metric = metric - min(metric);%transpose
                xs = metric(1:end-1).*cos(degBins(1:end-1));                % average
                ys = metric(1:end-1).*sin(degBins(1:end-1));               
                
                ang_hd = atan2(mean(ys),mean(xs));                          % mean direction
                
                mr = (cos(ang_hd)*sum(xs) + sin(ang_hd)*sum(ys)) / sum(metric(1:end-1)); % mean resultant length
                
                ind_fire = NaN(length(firing(:,1)),length(firing(1,:)));    %Indices of neuron activity/firing            
                MRLsave((length(firing(1,:))*itteration)-(length(firing(1,:))-cellNum),1) = mr; %save the MRL 
            end
        end
    end
end
out.mrall = MRLsave;
out.percentil99th = prctile(MRLsave,99);
out.firing = firing;
out.frequency = freq;
save('Shuffling_MLSpikeFULL.mat','out');
figure
histogram(MRLsave)
title('Shuffled MRL distribution')
ylabel('Number of Cells')
xlabel('Mean Resultant Length')
vline(out.percentil99th);
savefig('Shuffled_MRL_Hist_MLSpikeFULL.fig')
end

%% Subfunctions

%This function calculates the distance from the animal to boundaries of the environment at each behavioral data point.
%The distance calculation has to be done for all orientations around the animal centered on the animal’s
%current heading direction. That is to say that the animal’s current heading is always 0° and the distance
%to the boundaries is calculated for each of the 360 one-degree bins around the animal.
function [dis, ex, ey] = subfunc(rx,ry,hd, QP, degSamp)

mxd = sqrt((max(rx)-min(rx))^2 + (max(ry)-min(ry))^2);                      %sets bin radial maximum
degs = deg2rad(-180:degSamp:180);                                           %creating degree bins in radians
hd = deg2rad(hd);                                                           %converting Head direction in radians

edg = splitter(QP);                                                         %split 
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

%%output
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