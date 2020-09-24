function [out, detailed]= EBC_Social_Shuffle(normNeuC,Angle,DistanceH,deconvolve,DistBinSize,DegBinSize,name)
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

normNeuC = normNeuC';
FOVsize = 60;


%% Get behavior information
degBins = (-180:DegBinSize:179);                                  %Angle bins for EBC metric polar plot
degBins = degBins';                                     %reorient bins
degBins = deg2rad(degBins);                             %Convert to radians
mrall = zeros(1,length(normNeuC(1,:)));             %Firing frequency for each cell;
distanceBins = 0:DistBinSize:FOVsize;
fps = 30;                                               %Frames per second
spf = 1/fps;                                            %Seconds per frame
minDist = [1:2];
cutoff = 0;
timevar = [];

set(groot, 'defaultFigureUnits','normalized');
set(groot, 'defaultFigurePosition',[0 0 1 1]);

%% Calculate distances
[dis] = subfunc(DistanceH,Angle,1);   %calls funtion to bring back wall distances when neuron fired

firing = normNeuC;
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
                if i == 2 && i == mins
                    fire2 = fire((1800*(i-1)+1:end),1);
                    iTime2 = (1800*(i-1)+1:length(fire));
                    dis2 = dis((1800*(i-1)+1:end),:);
                elseif i ==2
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
                if i == 1 && i == mins
                    fire1 = fire(1:end,1);
                    iTime1 = (1:length(fire));
                    dis1 = dis(1:end,:);
                elseif i == 1
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
        if mins > 1
            ifire2 = find(fire2);
        else
            ifire2 = [];
        end                                     %Find indices for all non-zero values
        
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
                
                occ = occ(:,1:end-1); occ=occ';
                cutout = find(occ<cutoff);
                occ(cutout) = 0;
                nspk = nspk(:,1:end-1); nspk=nspk';
                %first half
                occ1 = occ1(:,1:end-1); occ1=occ1';
                %             cutout1 = find(occ<50);
                occ1(cutout) = 0;
                nspk1 = nspk1(:,1:end-1); nspk1=nspk1';
                %second half
                occ2 = occ2(:,1:end-1); occ2=occ2';
                %             cutout2 = find(occ<50);
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
                %             rm = (nspk./occ) * fps;
                nd = numel(thetaBins3d);
                rm = [rm rm rm];
                rm = CMBHOME.Utils.SmoothMat(rm,smooth(1:3),smooth(4));   % Smooth it
                rm = rm(:,nd+1:2*nd); % bring it back
                %first half
                rm1 = [rm1 rm1 rm1];
                rm1 = CMBHOME.Utils.SmoothMat(rm1,smooth(1:3),smooth(4));   % Smooth it
                rm1 = rm1(:,nd+1:2*nd); % bring it back
                %second half
                rm2 = [rm2 rm2 rm2];
                rm2 = CMBHOME.Utils.SmoothMat(rm2,smooth(1:3),smooth(4));   % Smooth it
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

% if find(MRLsave>0.6)
%     beep
%     pause
% end
out.mrall = MRLsave(:);
out.mrall1 = MRLsave1(:);
out.mrall2 = MRLsave2(:);
out.percentil99th = prctile(out.mrall,99);
out.firing = firing;
out.correlationEO = CORRsave(:);

save(name,'out');
figure
histogram(out.mrall)
title('Social MRL shuffle Distribution')
ylabel('Number of Cells')
xlabel('Mean Resultant Length')
savefig('Shuffled_MRL_Hist_Circ_CorrDonutSocial.fig')
end

%% Subfunctions

%This function calculates the distance from the animal to boundaries of the environment at each behavioral data point.
%The distance calculation has to be done for all orientations around the animal centered on the animal’s
%current heading direction. That is to say that the animal’s current heading is always 0° and the distance
%to the boundaries is calculated for each of the 360 one-degree bins around the animal.
function [dis] = subfunc(distance,hd, degSamp)
hd = hd + 180;
degs = (-180:degSamp:180);
dis = zeros(numel(distance),360);
dir = dis;

for i = 1:length(distance)
    dis(i,round(hd(i))) = distance(i);
end

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