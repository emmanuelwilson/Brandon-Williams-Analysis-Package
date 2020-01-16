function [out, detailed]= OccupancyWithBar(HD,tracking, pixX,ebcstats, varargin)
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

FOVsize = 21;

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
distanceBins = 0:1:FOVsize;                                  %set to look at half the length of the field which in our case is ~38cm (37.5 rounded up)
counter = 0 ;                                           %counter
fps = 30;                                               %Frames per second

%% Get structure of environment
%Identify where the bounds of the environment are located. Through a subfunction that allows the user to click
%on a plot of the positional data to indicate where the corners of the environment are located.

%% Calculate distances
if isempty(ebcstats)
    QP = findEdges(tracking);
else
    QP = ebcstats.QP;
end
degSamp = 1;                                                            %angle resolution
[dis, ex, ey] = subfunc(tracking(:,1),tracking(:,2),HD, QP, degSamp);   %calls funtion to bring back wall distances when neuron fired
dis_raw = dis;
dis = fillmissing(dis,'pchip',2);                                          %interpolates missing values
dis = dis*pixX;                                                             %Converts boundary distances from pixels to cm.
dis = circshift(dis,90,2);                                                  %shifts values by 90 degrees

%% Calculate raw maps:
thetaBins = deg2rad(linspace(-180,180,size(dis,2)));                    %angle bins
occ = NaN(length(thetaBins), length(distanceBins));                     %wall occupancy bins
distanceBins(end+1) = Inf;                                              %Adds an Infinity value at the end of the bins as safety procaution/break point
for i = 1:length(thetaBins)
    t = dis(:,i); %boundary distance for a particular bin
    for k = 1:length(distanceBins)-1
        inds = t>=distanceBins(k) & t<distanceBins(k+1);                %filter through the boundary distances
        occ(i,k) = sum (inds);                                          %Wall occupancy definition
    end
end
distanceBins = distanceBins(1:end-1);                                   %itteration through bins

% bring back to original dims
occ = occ(:,1:end-1); occ=occ';
cutout = find(occ<cutoff);
occ(cutout) = 0;

%% Smoothing
occ = [occ occ occ];
nd = numel(thetaBins);
occ = CMBHOME.Utils.SmoothMat(occ, smooth(1:2), smooth(3));
occ = occ(:, nd+1:2*nd);
occ(cutout) = 0;
occ = fliplr(occ);
occ = occ/fps;

%% Plots
figure(1);
clf;
n=1;
c = 1;
% Occupancy square
subplot(n,2,c);c=c+1;
imagesc(thetaBins,distanceBins,occ);
set(gca,'YDir','Normal'); colormap(jet);
caxis([0 101])
freezeColors
title('Occ')

% occupancy circular
subplot(n,2,c);c=c+1;
% the +pi/2 brings "forwards" to "up"
[t2, r2] = meshgrid(wrapTo2Pi(thetaBins+pi/2), distanceBins(1:end-1));
[x, y] = pol2cart(t2,r2);
surface(x,y, occ), shading interp
hold on
set(gca,'XTick',[],'YTick',[])
axis square
axis off
colormap(jet)
colorbar
% caxis([0 101])
%             set(gca, 'YDir','Normal','CLim',[0 prctile(occ(:),99)])
set(gca,'YDir','Normal')
freezeColors
title('occ')

%Save Results
saveas(gcf,'Occupancy.jpg'); %saving figure as a picture file (.jpg) in the new folder "EBCresults"
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