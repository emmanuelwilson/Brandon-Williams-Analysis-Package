function [wShift,hShift] = msAlignBetweenSessions2020(sessions,combs)
    
    sat = true;
    n = round(length(sessions)/5);
    if n*5 > length(sessions)
        m = 5;
    else
        m = 6;
    end
    
    figure
    for i = 1 : length(sessions)
        frame = sessions{i}.meanFiltFrame;
        subplot(m,n,i)
        pcolor(frame)
        shading flat
        daspect([1 1 1])
        colormap gray
    end
    %landmark selection
    curserW = zeros(length(sessions),1);
    curserH = zeros(length(sessions),1);
    count = 1;
    while sat  
        refFrame = sessions{count}.meanFiltFrame;
%         hLarge = fspecial('average', 80);                                      %Large Averaging filter
%         hSmall = fspecial('average', 10);                                       %Small Averaging filter
        
        figure(103)
        pcolor(refFrame)
        shading flat
        daspect([1 1 1])
        colormap gray
        title('Select landmark')
        [curserW(count), curserH(count), curserB] = ginput(1);
        
%         prompt = 'Do you want to keep these results? Y/N/S: ';
%         str = input(prompt,'s');
%         if isempty(str) || str == 'Y' || str == 'y'            
            count = count + 1;
%         elseif str == 'S' || str == 's'
%             wShift = 0;%pwShift;
%             hShift = 0;%phShift;
%             sat = false;
%         end
        if count == length(sessions)+1
            sat = false;
        end
    end
    
    %ROI selection.
    wShift = zeros(length(sessions));
    hShift = zeros(length(sessions));
%     rsh = 5;
%     rsw = 5;
    for i = 1 : length(combs)
%         frame = sessions{combs(i,1)}.meanFiltFrame;
%         refFrame = sessions{combs(i,2)}.meanFiltFrame;
        dw = curserW(combs(i,2)) - curserW(combs(i,1));
        dh = curserH(combs(i,2)) - curserH(combs(i,1));
%         fsize = size(frame);                
        itter = 150;
        ru = 25;
        rl = 3;
        rect = cell(itter,1);
        rectref = cell(itter,1);
        randvals = rl + (ru-rl).*rand(itter,1);
%         tformM = cell(itter,1);
%         movingRegistered = tformM;
        t1 = 100;
        t2 = 100;
        for j = 1: itter
            frame = sessions{combs(i,2)}.meanFiltFrame;
            refFrame = sessions{combs(i,1)}.meanFiltFrame;
            fsize = size(frame);
            if fsize(2) > curserW(combs(i,2)) + round(fsize(2)/randvals(j))
                rect{j}(2) = curserW(combs(i,2)) + round(fsize(2)/randvals(j));
                adjdiffup = 0;
            else
                adjdiffup = curserW(combs(i,2)) + round(fsize(2)/randvals(j)) - fsize(1);
                rect{j}(2) = curserW(combs(i,2)) + round(fsize(2)/randvals(j)) - adjdiffup;
            end
            if 1 < curserW(combs(i,2)) - round(fsize(1)/randvals(j)) - adjdiffup
                rect{j}(1) = curserW(combs(i,2)) - round(fsize(2)/randvals(j)) - adjdiffup;
            else
                adjdiffdn = 1 - curserW(combs(i,2)) - round(fsize(2)/randvals(j)) - adjdiffup;
                rect{j}(1) = curserW(combs(i,2)) - round(fsize(2)/randvals(j)) - adjdiffup + adjdiffdn;
            end
            
            if fsize(1) > curserH(combs(i,2)) + round(fsize(1)/randvals(j))
                rect{j}(4) = curserH(combs(i,2)) + round(fsize(1)/randvals(j));
                adjdiffr = 0;
            else
                adjdiffr = curserH(combs(i,2)) + round(fsize(1)/randvals(j)) - fsize(1);
                rect{j}(4) = curserH(combs(i,2)) + round(fsize(1)/randvals(j)) - adjdiffr;
            end
            if 1 < curserH(combs(i,2)) - round(fsize(1)/randvals(j)) - adjdiffr
                rect{j}(3) = curserH(combs(i,2)) - round(fsize(1)/randvals(j)) - adjdiffr;
            else
                adjdiffl = 1 - curserH(combs(i,2)) - round(fsize(1)/randvals(j)) - adjdiffr;
                rect{j}(3) = curserH(combs(i,2)) - round(fsize(1)/randvals(j)) - adjdiffr + adjdiffl;
            end
            
            if fsize(2) > curserW(combs(i,1)) + round(fsize(2)/randvals(j))
                rectref{j}(2) = curserW(combs(i,1)) + round(fsize(2)/randvals(j));
                adjdiffup = 0;
            else
                adjdiffup = curserW(combs(i,1)) + round(fsize(2)/randvals(j)) - fsize(1);
                rectref{j}(2) = curserW(combs(i,1)) + round(fsize(2)/randvals(j)) - adjdiffup;
            end
            if 1 < curserW(combs(i,1)) - round(fsize(1)/randvals(j)) - adjdiffup
                rectref{j}(1) = curserW(combs(i,1)) - round(fsize(2)/randvals(j)) - adjdiffup;
            else
                adjdiffdn = 1 - curserW(combs(i,1)) - round(fsize(2)/randvals(j)) - adjdiffup;
                rectref{j}(1) = curserW(combs(i,1)) - round(fsize(2)/randvals(j)) - adjdiffup + adjdiffdn;
            end
            
            if fsize(1) > curserH(combs(i,1)) + round(fsize(1)/randvals(j))
                rectref{j}(4) = curserH(combs(i,1)) + round(fsize(1)/randvals(j));
                adjdiffr = 0;
            else
                adjdiffr = curserH(combs(i,1)) + round(fsize(1)/randvals(j)) - fsize(1);
                rectref{j}(4) = curserH(combs(i,1)) + round(fsize(1)/randvals(j)) - adjdiffr;
            end
            if 1 < curserH(combs(i,1)) - round(fsize(1)/randvals(j)) - adjdiffr
                rectref{j}(3) = curserH(combs(i,1)) - round(fsize(1)/randvals(j)) - adjdiffr;
            else
                adjdiffl = 1 - curserH(combs(i,1)) - round(fsize(1)/randvals(j)) - adjdiffr;
                rectref{j}(3) = curserH(combs(i,1)) - round(fsize(1)/randvals(j)) - adjdiffr + adjdiffl;
            end
            
            ROI = uint16(rect{j});
            
%             refRect = [rect(3:4,j)- dh ; rect(1:2,j) - dw];
            refRect = uint16(rectref{j});
            refROI = uint16(refRect);
            
%             refFrame = (filter2(hSmall,refFrame) - filter2(hLarge, refFrame));
            
            
            refFrame = refFrame(refROI(3):refROI(4),refROI(1):refROI(2));
            refFrame = (refFrame-min(min(refFrame)))/max(max(refFrame-min(min(refFrame))));
            
%             frame = (filter2(hSmall,frame) - filter2(hLarge, frame));
            frame = frame(ROI(3):ROI(4),ROI(1):ROI(2));
            frame = (frame-min(min(frame)))/max(max(frame-min(min(frame))));
            if (length(frame(:,1)) < 16) || (length(frame(1,:)) < 16)                
                frame = padarray(frame,[7 7],0,'both');
                refFrame = padarray(refFrame,[7 7],0,'both');
            end
            ft = frame;
            refft = refFrame;
            
            [optimizer,metric] = imregconfig('multimodal');
            optimizer.InitialRadius = 0.001;
%             optimizer.Epsilon = 1.5e-4;
%             optimizer.GrowthFactor = 1.01;
%             optimizer.MaximumIterations = 500;
%             metric.NumberOfSpatialSamples = 1000;
%             metric.NumberOfHistogramBins = 50;            
            tformMtemp = imregtform(frame,refFrame,'translation',optimizer,metric); %rigid, similarity
            movingRegistered = imwarp(frame,tformMtemp,'OutputView',imref2d(size(refFrame)));
            tformMtemp.T = tformMtemp.T - [0,0,0;0,0,0;dw,dh,0];
            
            t1temp = tformMtemp.T(3,1);
            t2temp = tformMtemp.T(3,2);            
            
            if (abs(t1temp) < abs(t1) && abs(t2temp) < abs(t2)) || (abs(t1temp) < abs(t1) && abs(t2temp) == abs(t2)) || (abs(t1temp) == abs(t1) && abs(t2temp) < abs(t2))
                t1 = t1temp;
                t2 = t2temp;                
                refFrame = refft;
                tformM = tformMtemp;
                frame = ft;
                refFrame = refft;
                movereg = movingRegistered;
            end
        end
        
%         tformM.T(3,1) = -tformM.T(3,1);
%         tformM.T(3,2) = -tformM.T(3,2);
        
        subplot(2,3,1)
        imshow(refFrame);
        subplot(2,3,4)
        imshow(movereg);
        subplot(2,3,2)
        imshow(uint8(sessions{combs(i,1)}.meanFiltFrame));
        hold on
        plot(curserW(1),curserH(1),'+r','markersize',40);
        hold off
        subplot(2,3,5)
        imshow(uint8(sessions{combs(i,2)}.meanFiltFrame));
        hold on
        plot(curserW(1)-tformM.T(3,1),curserH(1)-tformM.T(3,2),'+r','markersize',40);
        hold off
        
        seg = imwarp(sessions{combs(i,2)}.outlines,tformM,'OutputView',imref2d(size(sessions{combs(i,1)}.meanFiltFrame)));
        
        subplot(2,3,[ 3  6])
        imshowpair(seg,sessions{combs(i,1)}.outlines)
        title(['wShift: ' num2str(tformM.T(3,1)) ' | hShift: ' num2str(tformM.T(3,2))]);
%         if isfield(vidObjFixed,'dateNum')
%             vidObjMoving.alignedDateNum = vidObjFixed.dateNum;
%         end
        wShift(combs(i,1),combs(i,2)) = tformM.T(3,1);% + pwShift;
        hShift(combs(i,1),combs(i,2)) = tformM.T(3,2);% + phShift;        
    end
    %}
end

