function [out, detailed]= msEgoCentricRateMapSplitEvenOddParallel_Social(normNeuC,Angle,DistanceH,deconvolve, cellIndex,DistBinSize,DegBinSize,name)
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
% name = ['EBCresults_Social',num2str(round(DistBinSize)),'A',num2str(DegBinSize),'smooth4_1_3_notSpike_notnorm'];
%% Get behavior information
ratemaps = zeros(FOVsize/DistBinSize,360/DegBinSize,length(normNeuC(1,:)));        %Probability ratemap values
ratemaps1 = zeros(FOVsize/DistBinSize,360/DegBinSize,length(normNeuC(1,:)));        %Probability ratemap values
ratemaps2 = zeros(FOVsize/DistBinSize,360/DegBinSize,length(normNeuC(1,:)));        %Probability ratemap values

degBins = (-180:DegBinSize:179);                                  %Angle bins for EBC metric polar plot
degBins = degBins';                                     %reorient bins
degBins = deg2rad(degBins);                             %Convert to radians
mrall = zeros(1,length(normNeuC(1,:)));             %Firing frequency for each cell;                                       %MRL for each cell
distanceBins = 0:DistBinSize:FOVsize;                                  %set to look at half the length of the field which in our case is ~38cm (37.5 rounded up)
mkdir(name)                          %Create new folder within current directory
counter = 0 ;                                           %counter
counter1 = 0;
counter2 = 0;
fps = 30;                                               %Frames per second
spf = 1/fps;                                            %Seconds per frame
minDist = [1:2];
cutoff = 0;
timevar = [];


set(groot, 'defaultFigureUnits','normalized');
set(groot, 'defaultFigurePosition',[0 0 1 1]);
%% Calculate distances
[dis] = subfunc(DistanceH,Angle,1);   %calls funtion to bring back wall distances when neuron fired

processed = false;
%Loop through every cell, extract and analize firing instances and boundary locations
for ctrack = 1 : length(cellIndex)
    cellNum = cellIndex(ctrack);        
    stime = tic;
    distanceBinsPar = distanceBins;
    if deconvolve == 1
        firing = normNeuC(:,cellNum);       
    else
        firing = normNeuC(:,cellNum);
    end
    fire = firing;                                  %Duplicate
    fire = fire - min(fire);
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
    ifire1 = find(fire1);
    if mins > 1
        ifire2 = find(fire2);
    else
        ifire2 = [];
    end
    if(~isempty(ifire) && ~isempty(ifire1) && ~isempty(ifire2))
        %full run
        ind_fire = ifire;                      %Add firing index to ms struct
        HDfiring = Angle(ifire);        %Head Direction of mouse at time of neural firing
        for j = 1 : length(ifire)
            %first half
            if j < length(ifire1)
                ind_fire1(j) = ifire1(j);                      %Add firing index to ms struct                
                HDfiring1(j) = Angle(iTime1(ifire1(j)));        %Head Direction of mouse at time of neural firing
            end
            %second half
            if j < length(ifire2)
                ind_fire2(j) = ifire2(j);                      %Add firing index to ms struct                
                HDfiring2(j) = Angle(iTime2(ifire2(j)));        %Head Direction of mouse at time of neural firing
            end
        end     
        
        %% Calculate raw maps:
        thetaBins = deg2rad(linspace(-180,180,size(dis,2)));                    %angle bins
        thetaBins3d = deg2rad(linspace(-180,180,round(size(dis,2)/DegBinSize)));                    %angle bins
        occ = NaN(length(thetaBins), length(distanceBinsPar));                     %wall occupancy bins
        nspk = occ;                                                             %Number of spikes bins
        nspk1 = nspk;
        nspk2 = nspk;
        occ1 = occ;
        occ2 = occ;
        distanceBinsPar(end+1) = Inf;                                              %Adds an Infinity value at the end of the bins as safety procaution/break point
%         distanceBins = distanceBins;
%         distanceBins = distanceBins;
        ci = ind_fire(:);                                            %firing instances of the cell
        ci1 = ind_fire1(:);                                            %firing instances of the cell
        ci2 = ind_fire2(:);                                            %firing instances of the cell        
        for i = 1:length(thetaBins)
            t = dis(:,i); %boundary distance for a particular bin
            t1 = dis1(:,i);
            t2 = dis2(:,i);
            for k = 1:length(distanceBinsPar)-1
                inds = t>=distanceBinsPar(k) & t<distanceBinsPar(k+1);                %filter through the boundary distances
                occ(i,k) = sum(inds);                                          %Wall occupancy definition
                inds = find(inds);                                              %find all non-zero boundary distances indices
                nspk(i,k) = sum(fire(intersect(inds,ci)));                      %Number of spike instances definition
                %first half
                inds1 = t1>=distanceBinsPar(k) & t1<distanceBinsPar(k+1);                %filter through the boundary distances
                occ1(i,k) = sum(inds1);                                          %Wall occupancy definition
                inds1 = find(inds1);                                              %find all non-zero boundary distances indices
                nspk1(i,k) = sum(fire1(intersect(inds1,ci1)));                      %Number of spike instances definition
                %second half
                inds2 = t2>=distanceBinsPar(k) & t2<distanceBinsPar(k+1);                %filter through the boundary distances
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
        
        distanceBinsPar = distanceBinsPar(1:end-1);                                   %itteration through bins
        if any(nspk(:)>0) && any(nspk1(:)>0) && any(nspk2(:)>0)
            counter = counter + 1;
            counter1 = counter1 +1;
            counter2 = counter2 +1;
            processed = true;
            
            % bring back to original dims
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
            maxrm1 = max(max(rm1));
            maxrm2 = max(max(rm2));
            if maxrm1>maxrm2
                maxrm = maxrm1;
            else
                maxrm = maxrm2;
            end
            
            %% Plots
            %The first three plots are rectangular versions of the three color-coded plots
            %Figure 7 is part of our work in progress for estimating a cell’s preferred distance.
            %Finally, figure 8 is a trajectory plot with heading color coded spike dots
            
            figure;
            set(0,'DefaultFigureWindowStyle' , 'normal')
            figure('Position', get(0, 'Screensize'));  
            figure('visible','off');
            n=3;
            m = 3;
            c = 1;

            % ratemap circular
            subplot(n,m,c);c=c+1;
            % the +pi/2 brings "forwards" to "up"
            [t2, r2] = meshgrid(wrapTo2Pi(thetaBins3d+pi/2), distanceBinsPar(1:end-1));
            [x, y] = pol2cart(t2,r2);
            h=surface(x,y, rm); shading interp
            
            hold on
            set(gca,'XTick',[],'YTick',[])
            axis square
            axis off
            colormap(jet)
%             set(gca, 'YDir','Normal','CLim',[0 prctile(rm(:), 99)])
            set(gca,'YDir','Normal')
            freezeColors
            title('rm')
            info = ['Max: ' num2str(max(max(rm)))];
            text(30,8,info);
            
%             %Trajectory map
%             if deconvolve == 1
%                 subplot(n,m,c);c=c+1;
%                 edg = splitter(QP);
%                 hold on
%                 plot(pixX*tracking(:,1),-pixY*tracking(:,2),'Color',[.7 .7 .7])
%                 colormap(hsv)
%                 xlim(pixX*[min(tracking(:,1)) max(tracking(:,1))]);ylim(pixY*[-max(tracking(:,2)) -min(tracking(:,2))])
%                 cx=pixX*cell_x(:);
%                 cy=pixY*cell_y(:);
%                 axis off
%                 axis square
%                 scatter(cx,-cy,38,HDfiring(:),'filled')
%                 set(gca,'YDir','Normal')
%                 caxis([0 360])
%                 title('Traj')                
%             else
%                 subplot(n,m,c);c=c+1;
%                 edg = splitter(QP);
%                 hold on
%                 plot(pixX*tracking(:,1),-pixY*tracking(:,2),'Color',[.7 .7 .7])
%                 colormap(hsv) 
%                 xlim(pixX*[min(tracking(:,1)) max(tracking(:,1))]);ylim(pixY*[-max(tracking(:,2)) -min(tracking(:,2))])                               
%             end
%             set(gca,'YDir','Normal')    
%             caxis([0 360])
%             title('Traj')
%             axis off
%             axis square
            
%             %EBC Metric            
            subplot(n,m,c);c=c+1;
            metric = mean(rm,1)';            
            if ~deconvolve                
                metric = metric - min(metric);
            end
            polarplot(degBins,metric)
            
            xs = metric(1:end-1).*cos(degBins(1:end-1)); % average
            ys = metric(1:end-1).*sin(degBins(1:end-1));
            
            coordlims=axis;
            
            ang_hd = atan2(mean(ys),mean(xs)); % mean direction
            
            mr = (cos(ang_hd)*sum(xs) + sin(ang_hd)*sum(ys)) / sum(metric(1:end-1)); % mean resultant length
            
            mag_hd = sqrt(sum(ys)^2+sum(xs)^2)/sqrt(sum(abs(ys))^2+sum(abs(xs))^2)*6.28; % for visualizations sake
            
            hold on;
            polarplot([ang_hd ang_hd ],[0 mr], 'r')
            pol = gca;
            pol.ThetaZeroLocation = 'top';
            pol.ThetaColor = 'none';
            pol.RTickLabel = [];
            hold off
            title('Wall Directionality')
            stat = ['MRL: ' num2str(mr) ' Angle : ' num2str(rad2deg(ang_hd))];
            text(0.2,coordlims(4),stat);              
            
            %FIRST HALF            
            % ratemap circular
            subplot(n,m,c);c=c+1;
            % the +pi/2 brings "forwards" to "up"
            [t2, r2] = meshgrid(wrapTo2Pi(thetaBins3d+pi/2), distanceBinsPar(1:end-1));
            [x, y] = pol2cart(t2,r2);
            h=surface(x,y, rm1); shading interp         
            
            hold on
            set(gca,'XTick',[],'YTick',[])
            axis square
            axis off
            colormap(jet)
            caxis([0 maxrm])
%             set(gca, 'YDir','Normal','CLim',[0 prctile(rm(:), 99)])
            set(gca,'YDir','Normal')
            freezeColors
            title('rm Odd Minutes')       
            info1 = ['Max: ' num2str(maxrm1) ' Correlation: ' num2str(corrpar)];
            text(30,coordlims(4),info1);
         
%             %Trajectory map
%             if deconvolve == 1
%                 subplot(n,m,c);c=c+1;
%                 edg = splitter(QP);
%                 hold on
%                 plot(pixX*tracking(frameMap(iTime1),1),-pixY*tracking(frameMap(iTime1),2),'Color',[.7 .7 .7])
%                 colormap(hsv)
%                 xlim(pixX*[min(tracking(:,1)) max(tracking(:,1))]);ylim(pixY*[-max(tracking(:,2)) -min(tracking(:,2))])
%                 cx=pixX*cell_x1(:);
%                 cy=pixY*cell_y1(:);
%                 axis off
%                 axis square
%                 scatter(cx,-cy,38,HDfiring1(:),'filled')
%                 set(gca,'YDir','Normal')
%                 caxis([0 360])
%                 title('Traj')          
%             else
%                 subplot(n,m,c);c=c+1;
%                 edg = splitter(QP);
%                 hold on
%                 plot(pixX*tracking(frameMap(ifire1),1),-pixY*tracking(frameMap(ifire1),2),'Color',[.7 .7 .7])
%                 colormap(hsv)
%                 xlim(pixX*[min(tracking(:,1)) max(tracking(:,1))]);ylim(pixY*[-max(tracking(:,2)) -min(tracking(:,2))])
%                 set(gca,'YDir','Normal')
%                 title('Traj')                   
%             end
%             set(gca,'YDir','Normal')    
%             caxis([0 360])
%             title('Traj')
%             axis off
%             axis square
            
%             %EBC Metric
            subplot(n,m,c);c=c+1;
            metric1 = mean(rm1,1)';            
            if ~deconvolve                
                metric1 = metric1 - min(metric1);
            end
            polarplot(degBins,metric1)
            
            xs = metric1(1:end-1).*cos(degBins(1:end-1)); % average
            ys = metric1(1:end-1).*sin(degBins(1:end-1));
            
            coordlims=axis;
            
            ang_hd1 = atan2(mean(ys),mean(xs)); % mean direction
            
            mr1 = (cos(ang_hd1)*sum(xs) + sin(ang_hd1)*sum(ys)) / sum(metric1(1:end-1)); % mean resultant length
            
            mag_hd1 = sqrt(sum(ys)^2+sum(xs)^2)/sqrt(sum(abs(ys))^2+sum(abs(xs))^2)*6.28; % for visualizations sake
            
            hold on;
            polarplot([ang_hd1 ang_hd1 ],[0 mr1], 'r')
            pol = gca;
            pol.ThetaZeroLocation = 'top';
            pol.ThetaColor = 'none';
            pol.RTickLabel = [];
            hold off
%             set(gca,'Color','none')
            title('Wall Directionality Obj: Odd Minutes')
            stat = ['MRL: ' num2str(mr1) 'Angle : ' num2str(rad2deg(ang_hd1))];
            text(0.2,coordlims(4),stat);
            
            %SECOND HALF           
            % ratemap circular
            subplot(n,m,c);c=c+1;
            % the +pi/2 brings "forwards" to "up"
            [t2, r2] = meshgrid(wrapTo2Pi(thetaBins3d+pi/2), distanceBinsPar(1:end-1));
            [x, y] = pol2cart(t2,r2);
            h=surface(x,y, rm2); shading interp

            hold on
            set(gca,'XTick',[],'YTick',[])
            axis square
            axis off
            colormap(jet)
            caxis([0 maxrm])
%             set(gca, 'YDir','Normal','CLim',[0 prctile(rm(:), 99)])
            set(gca,'YDir','Normal')
            freezeColors
            title('rm Even Minutes')
            info2 = ['Max: ' num2str(maxrm2)];
            text(30,coordlims(4),info2);
            
%             %Trajectory map
%             if deconvolve == 1
%                 subplot(n,m,c);c=c+1;
%                 edg = splitter(QP);
%                 hold on
%                 plot(pixX*tracking(frameMap(iTime2),1),-pixY*tracking(frameMap(iTime2),2),'Color',[.7 .7 .7])
%                 colormap(hsv)
%                 xlim(pixX*[min(tracking(:,1)) max(tracking(:,1))]);ylim(pixY*[-max(tracking(:,2)) -min(tracking(:,2))])
%                 cx=pixX*cell_x2(:);
%                 cy=pixY*cell_y2(:);
%                 axis off
%                 axis square
%                 scatter(cx,-cy,38,HDfiring2(:),'filled')
%                 set(gca,'YDir','Normal')
%                 caxis([0 360])
%                 title('Traj')                
%             else
%                 subplot(n,m,c);c=c+1;
%                 edg = splitter(QP);
%                 hold on
%                 plot(pixX*tracking(frameMap(ifire2),1),-pixY*tracking(frameMap(ifire2),2),'Color',[.7 .7 .7])
%                 colormap(hsv)
%                 xlim(pixX*[min(tracking(:,1)) max(tracking(:,1))]);ylim(pixY*[-max(tracking(:,2)) -min(tracking(:,2))])
%                 set(gca,'YDir','Normal')
%                 title('Traj')                 
%             end
%             set(gca,'YDir','Normal')    
%             caxis([0 360])
%             title('Traj')
%             axis off
%             axis square
            
%             %EBC Metric
            subplot(n,m,c);
            metric2 = mean(rm2,1)';
            if ~deconvolve                
                metric2 = metric2 - min(metric2);
            end
            polarplot(degBins,metric2)
            
            xs = metric2(1:end-1).*cos(degBins(1:end-1)); % average
            ys = metric2(1:end-1).*sin(degBins(1:end-1));
            
            coordlims=axis;
            
            ang_hd2 = atan2(mean(ys),mean(xs)); % mean direction
            
            mr2 = (cos(ang_hd2)*sum(xs) + sin(ang_hd2)*sum(ys)) / sum(metric2(1:end-1)); % mean resultant length
            
            mag_hd2 = sqrt(sum(ys)^2+sum(xs)^2)/sqrt(sum(abs(ys))^2+sum(abs(xs))^2)*6.28; % for visualizations sake
            
            hold on;
            polarplot([ang_hd2 ang_hd2 ],[0 mr2], 'r')
            pol = gca;
            pol.ThetaZeroLocation = 'top';
            pol.ThetaColor = 'none';
            pol.RTickLabel = [];
            hold off
            title('Wall Directionality: Even Minutes')
            stat = ['MRL: ' num2str(mr2) 'Angle : ' num2str(rad2deg(ang_hd2))];
            text(0.2,coordlims(4),stat);
                        
%             county  = county+1;
                       
            %Save Results
            saveas(gcf,[name,'/',num2str(cellNum),'EBC.jpg']); %saving figure as a picture file (.jpg) in the new folder "EBCresults"
            
            clf
            ind_fire = NaN(1,length(ind_fire(1,:))); %Indices of neuron activity/firing            
            HDfiring = ind_fire;                              %Indices of neuron activity for head direction 
            ind_fire1 = NaN(size(ind_fire1));            
%             cell_time1 = NaN(size(cell_time1));
            HDfiring1 = NaN;
            ind_fire2 = NaN(size(ind_fire2));            
%             cell_time2 = NaN(size(cell_time2));
            HDfiring2 = NaN;             
            
            close
        end
    end
    etime = toc(stime);
    timevar(ctrack) = etime;
    if processed
        ratemaps(:,:,ctrack) = rm;
        headangle(1,ctrack) = ang_hd;
        mrall(1,ctrack) = mr;
        
        ratemaps1(:,:,ctrack) = rm1;
        headangle1(1,ctrack) = ang_hd1;
        mrall1(1,ctrack) = mr1;
        
        ratemaps2(:,:,ctrack) = rm2;
        headangle2(1,ctrack) = ang_hd2;
        mrall2(1,ctrack) = mr2;
        
        correlations(ctrack) = corrpar;                
    end     
end
if processed
    out.time = timevar;
    out.corr = correlations;
    out.mravg = sum(mrall)/counter;
    out.mrall = mrall;
    out.mrall1 = mrall1;
    out.mrall2 = mrall2;
    out.headangle = headangle;
    out.headangle1 = headangle1;
    out.headangle2 = headangle2;
    out.rm = ratemaps;
    out.rm1 = ratemaps1;
    out.rm2 = ratemaps2;
    out.dis = dis;
    save([pwd,'/',name,'/','EBCstats.mat'],'out');
end
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