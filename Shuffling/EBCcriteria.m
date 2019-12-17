function [out, detailed]= EBCcriteria(ms,HD,tracking, frameMap, thresh, pixX, pixY, varargin)
%%Egocentric Boundary Cell criteria definition and filtering
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

%% Setup & Parse
p = inputParser;
p.addParameter('videoSamp', 1);                    % calculate every X frames of video
p.addParameter('degSamp', 1);                      % Degree bins
p.addParameter('heading', 1);                       % use heading (0--> Head direction)
p.addParameter('ifVideo', 0);                       % Play the video?
p.addParameter('labels', 1);                        % Plot forward?
p.addParameter('distanceBins', 0:1:38);           % How far to look (cm)
p.addParameter('boundaryMode', 1);                  % 0-> autolines, 1-> click, mat->useit
p.addParameter('ifLine', 0);
p.addParameter('figures',[0 0 0 1 1 1 1 1]);
p.addParameter('mergeFigures',1);                   % if 1 then puts in single figure
p.addParameter('smooth', [5 5 5])
p.parse(varargin{:});

%% Get behavior information

% ms.trace2 = ms.trace(1:length(frameMap)/2,1);
% ms.trace3 = ms.trace(length(frameMap)/2:length(frameMap),1);
% tracking2 = ms.trace2;
% tracking3 = tracking(length(tracking)/2:length(tracking),1:2);
% frameMap2 = frameMap(1:length(frameMap(:,1))/2,1);
% frameMap3 = frameMap(length(frameMap(:,1))/2:length(frameMap(:,1)),1);
% HD2 = HD(length(1:HD(:,1))/2,1);
% HD3 = HDdeg(length(HD(:,1))/2:length(HD(:,1)),1);
% ms.firing2 = ms.firing(1:length(ms.firing)/2,1);
% ms.firing3 = ms.firing(length(ms.firing)/2:length(ms.firing),1);

fps = 30;                                               %Frames per second
spf = 1/fps;                                            %Seconds per frame
ms.timestamp = frameMap.*spf;                           %time stamp in seconds
ms.ind_fire = NaN(ms.numFrames,length(ms.firing(1,:))); %Indices of neuron activity/firing
ms.cell_x = NaN(ms.numFrames,length(ms.firing(1,:)));   %X coordinate of locations where cell fired
ms.cell_y = NaN(ms.numFrames,length(ms.firing(1,:)));   %Y cooridnates of locations where cell fired
ms.cell_time = NaN(ms.numFrames,length(ms.firing(1,:)));%Time at when cell fired
ms.HDfiring = ms.ind_fire;                              %Indices of neuron activity for head direction
distanceBins = 0:1:37;                                  %set to look at half the length of the field which in our case is ~38cm (37.5 rounded up)
% mkdir EBCresultsTEST                                    %Create new folder within current directory
% direct = dir ('EBCresultsTEST');                         %Access new folder
degBins = (-90:2:269);                                  %Angle bins for EBC metric polar plot
degBins = degBins';
degBins = deg2rad(degBins);
mrtot = zeros(1,1,2);
freqFire = zeros(1,length(ms.firing(1,:)),2);
mrall = freqFire;
ratemaps = zeros(37,361,length(ms.firing(1,:)),2);
PrefDist = zeros(length(ms.firing(1,:)),1,2);
PrefAngl = PrefDist;
%% Get structure of environment
%Identify where the bounds of the environment are located. Through a subfunction that allows the user to click
%on a plot of the positional data to indicate where the corners of the environment are located.

QP = findEdges(tracking);       %finds the edge of the graph and mark the boundaries
tic
if length(frameMap) < length(HD)
    HD = HD(frameMap,:);   %Sink trace if not already sinked
    tracking = tracking(frameMap,:);
end

%% Half session seperation
firsthalf = length(HD(:,1))/2;
if firsthalf-length(frameMap)<0
    firsthalf = firsthalf +1;
    ms.firing = [zeros(1,length(ms.firing(1,:)));ms.firing];    
    tracking = [tracking(1,:); tracking];    
    HD = [HD(1,:);HD];
elseif firsthalf- length(frameMap) >0    
    ms.firing = [ms.firing;zeros(1,length(ms.firing(1,:)))];    
    tracking = [tracking;tracking(:,length(tracking(1,:)))];    
    HD = [HD;HD((length(HD(:,1))),:)];
end
firing(:,:,1) = ms.firing(1:1:firsthalf,:);
firing(:,:,2) = ms.firing(firsthalf:1:length(ms.firing),:);
HDeg(:,:,1) = HD(1:1:firsthalf);
HDeg(:,:,2) = HD(firsthalf:1:length(ms.firing));
track(:,:,1) = tracking(1:1:firsthalf,:);
track(:,:,2) = tracking(firsthalf:1:length(ms.firing),:);

%Loop through every cell, extract and analize firing instances and boundary locations
for itter = 1 :2
    %% Calculate distances
    degSamp = 1;                                                            %angle resolution
    [dis, ex, ey] = subfunc(track(:,1,itter),track(:,2,itter),HDeg(:,1,itter), QP, degSamp);   %calls funtion to bring back wall distances when neuron fired
    dis = fillmissing(dis,'spline',2);
    dis = dis*pixX;
    counter = 0;    
    for cellNum = 1 : length(ms.firing(1,:))
                
        fire = firing(:,:,itter);                                          %Firing trace for a half session
        
        fire(fire < thresh) = 0;                                           %apply threshold on firing
        ifire = find(fire);                                                %Find indices for all non-zero values
        if(~isempty(ifire))
            for j = 1 : length(ifire)
                ms.ind_fire(j,cellNum) = ifire(j);                          %Add firing index to ms struct
%                 ms.cell_x(j,cellNum) = tracking(ifire(j));                  %X postion of mouse during firing at sinked time
%                 ms.cell_y(j,cellNum) = tracking((ifire(j)),2);              %Y position of mouse during firing at sinked time
%                 ms.cell_time(j,cellNum) = ms.timestamp(ifire(j));           %Physiological time of firing
%                 ms.HDfiring(j,cellNum) = HDeg(ifire(j));                      %Head Direction of mouse at time of neural firing
            end            
            
            %% Calculate raw maps:
            thetaBins = deg2rad(linspace(-180,180,size(dis,2)));            %angle bins
            occ = NaN(length(thetaBins), length(distanceBins));             %wall occupancy bins
            nspk = occ;                                                     %Number of spikes bins
            distanceBins(end+1) = Inf;                                      %Adds an Infinity value at the end of the bins as safety procaution/break point
            ci = ms.ind_fire(:,cellNum);                                    %firing instances of the cell
            for i = 1:length(thetaBins)
                t = dis(:,i);                                               %boundary distance for a particular bin
                for k = 1:length(distanceBins)-1
                    inds = t>=distanceBins(k) & t<distanceBins(k+1);        %filter through the boundary distances
                    occ(i,k) = sum(inds);                                   %Wall occupancy definition
                    inds = find(inds);                                      %find all non-zero boundary distances indices
                    nspk(i,k) = length(intersect(inds,ci));                 %Number of spike instances definition
                end
            end
            distanceBins = distanceBins(1:end-1);                           %itteration through bins
            if any(nspk(:)>0)
                counter = counter + 1;
                
                % bring back to original dims
                occ = occ(:,1:end-1); occ=occ';
                nspk = nspk(:,1:end-1); nspk=nspk';
                
                %% Smoothing
                occ = [occ occ occ];
                nd = numel(thetaBins);
                occ = CMBHOME.Utils.SmoothMat(occ, smooth(1:2), smooth(3));
                occ = occ(:, nd+1:2*nd);
                
                nspk = [nspk nspk nspk];
                nspk = CMBHOME.Utils.SmoothMat(nspk,smooth(1:2),smooth(3)); % Smooth it
                nspk = nspk(:,nd+1:2*nd); % bring it back
                
                rm = (nspk./occ) * fps;
                rm = [rm rm rm];
                rm = CMBHOME.Utils.SmoothMat(rm,smooth(1:2),smooth(3));     % Smooth it
                rm = rm(:,nd+1:2*nd); % bring it back
                
                
                %% EBC Metric
                avgcount = zeros(1,i);
                metric = zeros(1,180);
                for it = 1 : i
                    avgcount(1,it) = mean(rm(:,it));
                    if mod(it,2) == 0
                        r = it/2;
                        metric(1,r) = (avgcount(1,it-1)+avgcount(1,it))/2;
                    end
                end
                metric = metric';
                
                xs = metric(1:end-1).*cos(degBins(1:end-1)); % average
                ys = metric(1:end-1).*sin(degBins(1:end-1));
                
                coordlims=axis;
                
                ang_hd = atan2(mean(ys),mean(xs)); % mean direction
                
                mr = (cos(ang_hd)*sum(xs) + sin(ang_hd)*sum(ys)) / sum(metric(1:end-1)); % mean resultant length
                
                mag_hd = sqrt(sum(ys)^2+sum(xs)^2)/sqrt(sum(abs(ys))^2+sum(abs(xs))^2)*6.28; % for visualizations sake
%                 hold on;
%                 polarplot([ang_hd ang_hd ],[0 mr], 'r')
%                 hold off
%                 title('Wall Directionality')
%                 stat = ['MRL: ' num2str(mr) 'Angle : ' num2str(rad2deg(ang_hd))];
%                 text(0.2,coordlims(4),stat);
                
                %Probability heatmap
                fig = figure;
                set(fig, 'Visible', 'off');
                [C,h] = contourf(rm,1);
                colormap gray
                set(gca,'Visible','off')                
%                 ax = gca;
%                 ax.Position = ax.OuterPosition;
                toplvl = zeros(size(rm));
                toplvl = (rm > h.LevelList(2));
                ydense = 0;
                ymax = 0;
                xdense = 0;
                xmax = 0;
                for yd = 1 : length(rm(1,:))
                    temp = sum(toplvl(:,yd));
                    if ymax < temp
                        ydense = yd;
                        ymax = temp;
                    end
                end
                for xd = 1 : length(rm(:,1))
                    temp = sum(toplvl(xd,:));
                    if xmax < temp
                        xdense = xd;
                        xmax = temp;
                    end
                end
                PrefDist(cellNum,1,itter) = ydense;
                PrefAngl(cellNum,1,itter) = xdense;
                
                %Save Results
%                 saveas(gcf,['EBCresultsTEST/',num2str(cellNum),'EBC.jpg']); %saving figure as a picture file (.jpg) in the new folder "EBCresults"
                ms.ind_fire = NaN(ms.numFrames,length(ms.firing(1,:))); %Indices of neuron activity/firing
%                 ms.cell_x = NaN(ms.numFrames,length(ms.firing(1,:)));   %X coordinate of locations where cell fired
%                 ms.cell_y = NaN(ms.numFrames,length(ms.firing(1,:)));   %Y cooridnates of locations where cell fired
%                 ms.cell_time = NaN(ms.numFrames,length(ms.firing(1,:)));%Time at when cell fired
%                 ms.HDfiring = ms.ind_fire;                              %Indices of neuron activity for head direction
                ratemaps(:,:,cellNum,itter) = rm;
                mrtot(1,1,itter) = mrtot(1,1,itter) + mr;
                freqFire(1,cellNum,itter) = length(ifire);
                mrall(1,cellNum,itter) = mr;               
                clf
            end
        end
    end
    out.mravg(1,itter) = mrtot(:,:,itter)/counter;
    out.mrall(1,:,itter) = mrall(:,:,itter);           
end
out.freqFire = freqFire;
out.rm = ratemaps;
out.PrefDist = PrefDist;
out.PrefAngl = PrefAngl;
toc
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