%%Will create a new folder and generate/save EBC plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This script is designed for 4 animals. Will identify which cells you want
%to keep and generate allocentric ratemap of activity, trajecotry plot with
%colourcoated points indicating acitivity location and orientation, MRL
%polar plot, Binarized EBC plot and Continous EBC plot. You need to specify
%the folder name in line 21 and line 511.

%INPUTS: These are variables which need to be in your workspace
% - indOut: index of cells to eliminate
% - ms: requires ms1-ms4 depending on how many animals you want to look at
% - mice: matrix containing the binarized EBC ratemap of all cells across
% animals
% - micec: same as mice but for continous signal
% -frameMap: requires frameMap1-frameMap4 depending on how many animals
% - SINKdata: requires SINKdata trajectory for each animal
% You will need access to the "Binarize.m" function 
%
%Output: Saved figure saved in the folder that you specified
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mkdir BESTebcs96thPer2deg
count = 1;
m = 2;
n = 2;

degBins = (-180:2:179);                                  %Angle bins for EBC metric polar plot
degBins = degBins';                                     %reorient bins
degBins = deg2rad(degBins);                             %Convert to radians

% px = 0;
% py = 0;
thetaBins = deg2rad(linspace(-180,180,361));
thetaBins2 = deg2rad(linspace(-180,180,181)); 
distanceBins = 0:1:30;
mice2 = zeros(length(mice(:,1,1)),181,length(mice(1,1,:)));
micec2 = zeros(length(mice(:,1,1)),181,length(mice(1,1,:)));

for i = 1: length(mice(:,1,1))
    for j = 1 : length(mice(1,:,1))
        if mod(j,2) == 0
            mice2(i,j/2,:) = mean(mice(i,j-1:j,:),2);
            micec2(i,j/2,:) = mean(micec(i,j-1:j,:),2);
        elseif j == length(mice(1,:,1))
            mice2(i,round(j/2),:) = mice(i,j,:);
            micec2(i,round(j/2),:) = micec(i,j,:);
        end
    end
end

firePeaks1 = Binarize(ms1);
firePeaks2 = Binarize(ms2);
firePeaks3 = Binarize(ms3);
firePeaks4 = Binarize(ms4);

probMap = zeros(30,30);                                 %Probablility open feild allocentric rate map
% J = length(probMap(:));
%%
for i = 612: length(mice(1,1,:))
    figure(1)
    probMap = zeros(30,30);
    if isempty(find(indOut == i))
        if i <= length(mouse1.mrall(1,:))                                   %First Animal
            pixX = 0.345;
            pixY = 0.35;
            firing = firePeaks1.binarizedTraces(:,i);%ms.FiltTraces(:,cellNum);%ms.firing(:,cellNum);                          %Extract firing trace
            firing = circshift(firing,-6,1);
            firing(end-6 : end, :) = 0;
            fire = firing;                                  %Duplicate
            if min(fire)<0
                fire = fire + abs(min(fire));
            else
                fire = fire - min(fire);
            end
            ifire = find(fire);
            
            %organizing mouse location
            lx =  max(mouse1.QP(:,1)) -  min(mouse1.QP(:,1));
            ly =  max(mouse1.QP(:,2)) -  min(mouse1.QP(:,2));
            locx = SINKdata1(frameMap1,1);
            locy = SINKdata1(frameMap1,2);
            locx = ((locx - min(mouse1.QP(:,1)))/lx)*30;
            locy = ((locy - min(mouse1.QP(:,2)))/lx)*30;
            difx = locx - round(locx);
            locx(difx>0) = round(locx(difx>0))+1;
            locx(difx<=0) = round(locx(difx<=0));
            dify = locy - round(locy);
            locy(dify>0) = round(locy(dify>0))+1;
            locy(dify<=0) = round(locy(dify<=0));
            
            for j = 1 : length(ifire)
                ind_fire(j,i) = ifire(j);                      %Add firing index to ms struct
                cell_x(j,i) = SINKdata1((frameMap1(ifire(j))),1);  %X postion of mouse during firing at sinked time
                cell_y(j,i) = SINKdata1((frameMap1(ifire(j))),2);%Y position of mouse during firing at sinked time
                %             cell_time(j,i) = ms1.timestamp(ifire(j));       %Physiological time of firing
                HDfiring(j,i) = HDdeg1(frameMap1(ifire(j)));        %Head Direction of mouse at time of neural firing
                
%                 if cell_x(j,i) < max(mouse1.QP(:,1)) && cell_x(j,i)> min(mouse1.QP(:,1))
%                     vx = ((cell_x(j,i)-min(mouse1.QP(:,1)))/lx)*30;
%                     px = vx - round(vx);
%                     if px>0
%                         px = round(vx) +1;
%                     else
%                         px = round(vx);
%                     end
%                 end
%                 if cell_y(j,i) < max(mouse1.QP(:,2)) && cell_y(j,i)> min(mouse1.QP(:,2))
%                     ly =  max(mouse1.QP(:,2)) -  min(mouse1.QP(:,2));
%                     vy = ((cell_y(j,i)-min(mouse1.QP(:,2)))/ly)*30;
%                     py = vy - round(vy);
%                     if py>0
%                         py = round(vy) +1;
%                     else
%                         py = round(vy);
%                     end
%                 end
%                 probMap(py,px) = probMap(py,px) +1;
            end
            
%             for j = 1 : length(probMap(:))
%                 [y, x] = ind2sub(size(probMap),j);
%                 if x == 16 && y == 14
%                     
%                 end
%                 if probMap(j) == 0 && isempty(find(locy(find(locx == x)) == y))
%                     probMap(j) = NaN;
%                     J(j) = j;
%                 end
%             end
%             J = unique(J);
%             if J(1) == 0
%                 J(1) = [];
%             end
%             probMap = CMBHOME.Utils.SmoothMat(probMap, smooth(1:2), smooth(3));
%             [nr, nc] = size(probMap);
%             subplot(m,n,count);
%             imAlpha=ones(size(probMap));
%             imAlpha(J)=0;
%             imagesc(probMap,'AlphaData',imAlpha);
%             set(gca,'color',0*[0 0 0]);
%             colormap(gca, jet)
%             title(['Allocentric RateMap' num2str(i)])
%             axis off
%             axis square
%             count = count+1;
            
            %Trajectory Plot
            edg = splitter(mouse1.QP);
            subplot(m,n,count);
            hold on
            plot(pixX*SINKdata1(:,1),-pixY*SINKdata1(:,2),'Color',[.7 .7 .7])
            colormap(gca,hsv)
            xlim(pixX*[min(SINKdata1(:,1)) max(SINKdata1(:,1))]);ylim(pixY*[-max(SINKdata1(:,2)) -min(SINKdata1(:,2))])
            cx=pixX*cell_x(:,i);
            cy=pixY*cell_y(:,i);
            scatter(cx,-cy,38,HDfiring(:,i),'filled')
            set(gca,'YDir','Normal')
            caxis([0 360])
            title(['Traj' num2str(i) ])
            axis off
            axis square
            hold off
            count = count+1;
            
        elseif i <= (length(mouse2.mrall(1,:)) + length(mouse1.mrall(1,:))) %Second Animal
            lx =  max(mouse2.QP(:,1)) -  min(mouse2.QP(:,1));
            ly =  max(mouse2.QP(:,2)) -  min(mouse2.QP(:,2));
            locx = SINKdata2(frameMap2,1);
            locy = SINKdata2(frameMap2,2);
            locx = ((locx - min(mouse2.QP(:,1)))/lx)*30;
            locy = ((locy - min(mouse2.QP(:,2)))/lx)*30;
            difx = locx - round(locx);
            locx(difx>0) = round(locx(difx>0))+1;
            locx(difx<=0) = round(locx(difx<=0));
            dify = locy - round(locy);
            locy(dify>0) = round(locy(dify>0))+1;
            locy(dify<=0) = round(locy(dify<=0));
            
            pixX = 0.345;
            pixY = 0.35;
            cellNum = i- length(mouse1.mrall(1,:));
            firing = firePeaks2.binarizedTraces(:,cellNum);%ms.FiltTraces(:,cellNum);%ms.firing(:,cellNum);                          %Extract firing trace
            firing = circshift(firing,-6,1);
            firing(end-6 : end, :) = 0;
            fire = firing;                                  %Duplicate
            if min(fire)<0
                fire = fire + abs(min(fire));
            else
                fire = fire - min(fire);
            end
            ifire = find(fire);
            for j = 1 : length(ifire)
                ind_fire(j,i) = ifire(j);                      %Add firing index to ms struct
                cell_x(j,i) = SINKdata2((frameMap2(ifire(j))),1);  %X postion of mouse during firing at sinked time
                cell_y(j,i) = SINKdata2((frameMap2(ifire(j))),2);%Y position of mouse during firing at sinked time
                %             cell_time(j,i) = ms1.timestamp(ifire(j));       %Physiological time of firing
                HDfiring(j,i) = HDdeg2(frameMap2(ifire(j)));        %Head Direction of mouse at time of neural firing
                
                if cell_x(j,i) < max(mouse2.QP(:,1)) && cell_x(j,i)> min(mouse2.QP(:,1))
                    lx =  max(mouse2.QP(:,1)) -  min(mouse2.QP(:,1));
                    vx = ((cell_x(j,i)-min(mouse2.QP(:,1)))/lx)*30;
                    px = vx - round(vx);
                    if px>0
                        px = round(vx) +1;
                    else
                        px = round(vx);
                    end
                end
                if cell_y(j,i) < max(mouse2.QP(:,2)) && cell_y(j,i)> min(mouse2.QP(:,2))
                    lx =  max(mouse2.QP(:,2)) -  min(mouse2.QP(:,2));
                    vy = ((cell_y(j,i)-min(mouse2.QP(:,2)))/lx)*30;
                    py = vy - round(vy);
                    if py>0
                        py = round(vy) +1;
                    else
                        py = round(vy);
                    end
                end
                probMap(py,px) = probMap(py,px) +1;
            end
%             for j = 1 : length(probMap(:))
%                 [y, x] = ind2sub(size(probMap),j);
%                 if probMap(j) == 0 && isempty(find(locy(find(locx == x)) == y))
%                     probMap(j) = NaN;
%                     J(j) = j;
%                 end
%             end
%             J = unique(J);
%             if J(1) == 0
%                 J(1) = [];
%             end
%             probMap = CMBHOME.Utils.SmoothMat(probMap, smooth(1:2), smooth(3));
%             [nr, nc] = size(probMap);
%             subplot(m,n,count);
%             imAlpha=ones(size(probMap));
%             imAlpha(J)=0;
%             imagesc(probMap,'AlphaData',imAlpha);
%             set(gca,'color',0*[0 0 0]);
%             colormap(gca, jet)
%             title(['Allocentric RateMap' num2str(i)])
%             axis off
%             axis square
%             count = count+1;
            %Trajectory Plot
            subplot(m,n,count);
            edg = splitter(mouse2.QP);
            hold on
            plot(pixX*SINKdata2(:,1),-pixY*SINKdata2(:,2),'Color',[.7 .7 .7])
            colormap(gca,hsv)
            xlim(pixX*[min(SINKdata2(:,1)) max(SINKdata2(:,1))]);ylim(pixY*[-max(SINKdata2(:,2)) -min(SINKdata2(:,2))])
            cx=pixX*cell_x(:,i);
            cy=pixY*cell_y(:,i);
            scatter(cx,-cy,38,HDfiring(:,i),'filled')
            set(gca,'YDir','Normal')
            caxis([0 360])
            title(['Traj' num2str(i) ])
            axis off
            axis square
            count = count+1;
            
        elseif i <= (length(mouse3.mrall(1,:)) + length(mouse2.mrall(1,:)) + length(mouse1.mrall(1,:)))
            lx =  max(mouse3.QP(:,1)) -  min(mouse3.QP(:,1));
            ly =  max(mouse3.QP(:,2)) -  min(mouse3.QP(:,2));
            frameMap3 = frameMap3(1:find(frameMap3 == length(SINKdata3(:,1))));
            locx = SINKdata3(frameMap3,1);
            locy = SINKdata3(frameMap3,2);
            locx = ((locx - min(mouse3.QP(:,1)))/lx)*30;
            locy = ((locy - min(mouse3.QP(:,2)))/lx)*30;
            difx = locx - round(locx);
            locx(difx>0) = round(locx(difx>0))+1;
            locx(difx<=0) = round(locx(difx<=0));
            dify = locy - round(locy);
            locy(dify>0) = round(locy(dify>0))+1;
            locy(dify<=0) = round(locy(dify<=0));
            cellNum = i- length(mouse1.mrall(1,:)) - length(mouse2.mrall(1,:));
            pixX = 0.345;
            pixY = 0.35;
            firing = firePeaks3.binarizedTraces(:,cellNum);%ms.FiltTraces(:,cellNum);%ms.firing(:,cellNum);                          %Extract firing trace
            firing = circshift(firing,-6,1);
            firing(end-6 : end, :) = 0;
            fire = firing;                                  %Duplicate
            if min(fire)<0
                fire = fire + abs(min(fire));
            else
                fire = fire - min(fire);
            end
            ifire = find(fire);
            for j = 1 : length(ifire)
                ind_fire(j,i) = ifire(j);                      %Add firing index to ms struct
                cell_x(j,i) = SINKdata3((frameMap3(ifire(j))),1);  %X postion of mouse during firing at sinked time
                cell_y(j,i) = SINKdata3((frameMap3(ifire(j))),2);%Y position of mouse during firing at sinked time
                %             cell_time(j,i) = ms1.timestamp(ifire(j));       %Physiological time of firing
                HDfiring(j,i) = HDdeg3(frameMap3(ifire(j)));        %Head Direction of mouse at time of neural firing
                
                if cell_x(j,i) < max(mouse3.QP(:,1)) && cell_x(j,i)> min(mouse3.QP(:,1))
                    lx =  max(mouse3.QP(:,1)) -  min(mouse3.QP(:,1));
                    vx = ((cell_x(j,i)-min(mouse3.QP(:,1)))/lx)*30;
                    px = vx - round(vx);
                    if px>0
                        px = round(vx) +1;
                    else
                        px = round(vx);
                    end
                end
                if cell_y(j,i) < max(mouse3.QP(:,2)) && cell_y(j,i)> min(mouse3.QP(:,2))
                    lx =  max(mouse3.QP(:,2)) -  min(mouse3.QP(:,2));
                    vy = ((cell_y(j,i)-min(mouse3.QP(:,2)))/lx)*30;
                    py = vy - round(vy);
                    if py>0
                        py = round(vy) +1;
                    else
                        py = round(vy);
                    end
                end
                probMap(py,px) = probMap(py,px) +1;
            end
            
%             for j = 1 : length(probMap(:))
%                 [y, x] = ind2sub(size(probMap),j);
%                 if probMap(j) == 0 && isempty(find(locy(find(locx == x)) == y))
%                     probMap(j) = NaN;
%                     J(j) = j;
%                 end
%             end
%             J = unique(J);
%             if J(1) == 0
%                 J(1) = [];
%             end
%             
%             probMap = CMBHOME.Utils.SmoothMat(probMap, smooth(1:2), smooth(3));
%             [nr, nc] = size(probMap);
%             subplot(m,n,count);
%             imAlpha=ones(size(probMap));
%             imAlpha(J)=0;
%             imagesc(probMap,'AlphaData',imAlpha);
%             set(gca,'color',0*[0 0 0]);
%             colormap(gca, jet)
%             title(['Allocentric RateMap' num2str(i)])
%             axis off
%             axis square
%             count = count+1;

            %Trajectory Plot
            subplot(m,n,count);
            edg = splitter(mouse3.QP);
            hold on
            plot(pixX*SINKdata3(:,1),-pixY*SINKdata3(:,2),'Color',[.7 .7 .7])
            colormap(gca,hsv)
            xlim(pixX*[min(SINKdata3(:,1)) max(SINKdata3(:,1))]);ylim(pixY*[-max(SINKdata3(:,2)) -min(SINKdata3(:,2))])
            cx=pixX*cell_x(:,i);
            cy=pixY*cell_y(:,i);
            scatter(cx,-cy,38,HDfiring(:,i),'filled')
            set(gca,'YDir','Normal')
            caxis([0 360])
            title(['Traj' num2str(i) ])
            axis off
            axis square
            count = count+1;
        elseif i <= (length(mouse4.mrall(1,:)) + length(mouse3.mrall(1,:)) + length(mouse2.mrall(1,:)) + length(mouse1.mrall(1,:)))
            lx =  max(mouse4.QP(:,1)) -  min(mouse4.QP(:,1));
            ly =  max(mouse4.QP(:,2)) -  min(mouse4.QP(:,2));
            locx = SINKdata4(frameMap4,1);
            locy = SINKdata4(frameMap4,2);
            locx = ((locx - min(mouse4.QP(:,1)))/lx)*30;
            locy = ((locy - min(mouse4.QP(:,2)))/lx)*30;
            difx = locx - round(locx);
            locx(difx>0) = round(locx(difx>0))+1;
            locx(difx<=0) = round(locx(difx<=0));
            dify = locy - round(locy);
            locy(dify>0) = round(locy(dify>0))+1;
            locy(dify<=0) = round(locy(dify<=0));
            cellNum = i- length(mouse1.mrall(1,:)) - length(mouse2.mrall(1,:)) - length(mouse3.mrall(1,:));
            pixX = 0.345;
            pixY = 0.35;
            firing = firePeaks4.binarizedTraces(:,cellNum);%ms.FiltTraces(:,cellNum);%ms.firing(:,cellNum);                          %Extract firing trace
            firing = circshift(firing,-6,1);
            firing(end-6 : end, :) = 0;
            fire = firing;                                  %Duplicate
            if min(fire)<0
                fire = fire + abs(min(fire));
            else
                fire = fire - min(fire);
            end
            ifire = find(fire);
            for j = 1 : length(ifire)
                ind_fire(j,i) = ifire(j);                      %Add firing index to ms struct
                cell_x(j,i) = SINKdata4((frameMap4(ifire(j))),1);  %X postion of mouse during firing at sinked time
                cell_y(j,i) = SINKdata4((frameMap4(ifire(j))),2);%Y position of mouse during firing at sinked time
                %             cell_time(j,i) = ms1.timestamp(ifire(j));       %Physiological time of firing
                HDfiring(j,i) = HDdeg4(frameMap4(ifire(j)));        %Head Direction of mouse at time of neural firing
                
%                 if cell_x(j,i) < max(mouse4.QP(:,1)) && cell_x(j,i)> min(mouse4.QP(:,1))
%                     lx =  max(mouse4.QP(:,1)) -  min(mouse4.QP(:,1));
%                     vx = ((cell_x(j,i)-min(mouse4.QP(:,1)))/lx)*30;
%                     px = vx - round(vx);
%                     if px>0
%                         px = round(vx) +1;
%                     else
%                         px = round(vx);
%                     end
%                 end
%                 if cell_y(j,i) < max(mouse4.QP(:,2)) && cell_y(j,i)> min(mouse4.QP(:,2))
%                     lx =  max(mouse4.QP(:,2)) -  min(mouse4.QP(:,2));
%                     vy = ((cell_y(j,i)-min(mouse4.QP(:,2)))/lx)*30;
%                     py = vy - round(vy);
%                     if py>0
%                         py = round(vy) +1;
%                     else
%                         py = round(vy);
%                     end
%                 end
%                 probMap(py,px) = probMap(py,px) +1;
            end
%             for j = 1 : length(probMap(:))
%                 [y, x] = ind2sub(size(probMap),j);
%                 if probMap(j) == 0 && isempty(find(locy(find(locx == x)) == y))
%                     probMap(j) = NaN;
%                     J(j) = j;
%                 end
%             end
%             J = unique(J);
%             if J(1) == 0
%                 J(1) = [];
%             end
%             probMap = CMBHOME.Utils.SmoothMat(probMap, smooth(1:2), smooth(3));
%             [nr, nc] = size(probMap);
%             subplot(m,n,count);
%             imAlpha=ones(size(probMap));
%             imAlpha(J)=0;
%             imagesc(probMap,'AlphaData',imAlpha);
%             set(gca,'color',0*[0 0 0]);
%             colormap(gca, jet)
%             title(['Allocentric RateMap' num2str(i)])
%             axis off
%             axis square
%             count = count+1;
            
            %Trajectory Plot
            subplot(m,n,count);
            edg = splitter(mouse4.QP);
            hold on
            plot(pixX*SINKdata4(:,1),-pixY*SINKdata4(:,2),'Color',[.7 .7 .7])
            colormap(gca,hsv)
            xlim(pixX*[min(SINKdata4(:,1)) max(SINKdata4(:,1))]);ylim(pixY*[-max(SINKdata4(:,2)) -min(SINKdata4(:,2))])
            cx=pixX*cell_x(:,i);
            cy=pixY*cell_y(:,i);
            scatter(cx,-cy,38,HDfiring(:,i),'filled')
            set(gca,'YDir','Normal')
            caxis([0 360])
            title(['Traj' num2str(i) ])
            axis off
            axis square
            count = count+1;
        else
            disp('ERROR: Unknown cell session location');
            
        end
        
        %EBC Metric and MRL polar plot        
        avgcount = zeros(1,length(thetaBins));
        metric = zeros(1,180);
        subplot(m,n,count);        
        %             contour(rm)
        r = 0;
        for it = 1 : length(thetaBins)
            avgcount(1,it) = mean(mice(:,it,i));
            if mod(it,2) == 0
                r = it/2;
                metric(1,r) = (avgcount(1,it-1)+avgcount(1,it))/2;
            end
        end
        metric = metric';
        
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
        hold off
        title('Wall Directionality')
        stat = ['MRL: ' num2str(mr) 'Angle : ' num2str(rad2deg(ang_hd))];
        text(0.2,coordlims(4),stat);        
        count = count +1;                     
        
        %Binarized EBC plot
        subplot(m,n,count)
        [t2, r2] = meshgrid(wrapTo2Pi(thetaBins2+pi/2), distanceBins(1:end-1));
        [x, y] = pol2cart(t2,r2);
        surface(x,y, mice2(:,:,i)); shading interp
        colormap(gca, jet)
        hold on
        set(gca,'XTick',[],'YTick',[])
        axis off
        axis square
        axis equal
        set(gca,'YDir','Normal')
        %         freezeColors
        title(['rm' num2str(i)])
        count = count +1;
        
        %Continuous EBC plot
        subplot(m,n,count)
        [t2, r2] = meshgrid(wrapTo2Pi(thetaBins2+pi/2), distanceBins(1:end-1));
        [x, y] = pol2cart(t2,r2);
        surface(x,y, micec2(:,:,i)); shading interp
        colormap(gca, jet)
        hold on
        set(gca,'XTick',[],'YTick',[])
        axis off
        axis square
        axis equal
        set(gca,'YDir','Normal')
        %         freezeColors
        title(['rm' num2str(i)])
        count = count +1;
    end
    locx = [];
    locy = [];
%     J = [];
    if count > 3
        %         print('FillPageFigure','-dpdf','-fillpage')
        saveas(gcf,['BESTebcs96thPer2deg/',num2str(i),'EBC.pdf']);
        clf
        count = 1;
    end
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