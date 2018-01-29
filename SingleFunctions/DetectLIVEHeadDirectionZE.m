clear;
clc;
close all;
objects = imaqfind; %find video input objects in memory

delete(objects) %delete a video input object from memory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of frames to inspect
FrameNum = 100;

REDco = zeros(FrameNum, 20); %Red Centroid coordinate storage
GREENco = zeros(FrameNum,20);%Green Centroid Coordinate storage
redThresh = 0.007;   % 0.01 * 0.0039
greenThresh = 0.007;  % 0.01 //best: 0.007
difference_threshold_red = 350;  % 250   //best: 350
difference_threshold_green = 700;  % 260  //best: 400  700
CentCo = [0,0];

vid = videoinput('winvideo',1,'YUY2_320x240');      % (Window video format, adaptor#, video resolution)
%%%PLEASE type in the resolution parameters below in order to enforce boundary conditions %%%
X_resolution = 320;
Y_resolution = 240;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%.
%%%%%
set(vid,'FramesPerTrigger',inf);
set(vid,'ReturnedColorspace','rgb')                 %Colour format
vid.FrameGrabInterval = 2;                         %interval of frame extraction ie. 5 = 1 frame every 5 milliseconds

start(vid);                                         %video aquiring process

while(vid.FramesAcquired <= FrameNum)                    %# of frames being processed
    
    frame = vid.FramesAcquired;
    ThisFrame=getsnapshot(vid); 
    
    diffFrameRed = imsubtract(thisFrame(:,:,1), rgb2gray(thisFrame)); % Get red component of the image
    diffFrameRed = medfilt2(diffFrameRed, [3 3]); % Filter out the noise by using median filter
    binFrameRed = imbinarize(diffFrameRed, redThresh); % Convert the image into binary image with the red objects as white
    
    CC_red = bwconncomp(binFrameRed);
    S_red = regionprops(CC_red, 'Centroid');
    
    centroids_red = cat(1, S_red.Centroid);
    
    if ~isempty(centroids_red)
        centroids_red_last_frame_previous_vid = centroids_red;
        if ind_vid == 1
            no_red_led = 0;
        end
    end 
    
    diffFrameGreen = imsubtract(thisFrame(:,:,2), rgb2gray(thisFrame)); % Get green component of the image
    diffFrameGreen = medfilt2(diffFrameGreen, [3 3]); % Filter out the noise by using median filter
    binFrameGreen = imbinarize(diffFrameGreen, greenThresh); % Convert the image into binary image with the green objects as white
    
    CC_green = bwconncomp(binFrameGreen);
    S_green = regionprops(CC_green, 'Centroid');
    
    centroids_green = cat(1, S_green.Centroid);
    
    if ~isempty(centroids_green)
        centroids_green_last_frame_previous_vid = centroids_green;
        if ind_vid == 1
            no_green_led = 0;
        end
    end
    
    if (frame == first_frame) && (no_red_led == 0) && (no_green_led ==0)
        
        imshow(thisFrame);
        hold on
        if (~isempty(position_track))
            plot(position_track(:,1), position_track(:,2), 'w')
        end
        hold on
        if (~isempty(centroids_red))
            plot(centroids_red(:,1), centroids_red(:,2), 'r*')
        else
            centroids_red = centroids_red_last_frame_previous_vid;
            plot(centroids_red(:,1), centroids_red(:,2), 'r*')
        end
        hold on
        if (~isempty(centroids_green))
            plot(centroids_green(:,1), centroids_green(:,2), 'g*')
        else
            centroids_green = centroids_green_last_frame_previous_vid;
            plot(centroids_green(:,1), centroids_green(:,2), 'g*')
        end
        %{
        if ind_vid == 1
            display('Choose red spot (use right-click)')
            [x_red, y_red] = getpts
            display('Choose green spot (use right-click)')
            [x_green, y_green] = getpts
            SquareDifferenceVect_red = [(centroids_red(:,1) - x_red).^2 + (centroids_red(:,2) - y_red).^2];
            SquareDifferenceVect_green = [(centroids_green(:,1) - x_green).^2 + (centroids_green(:,2) - y_green).^2];
            [valDiff_red idx_red_1] = min(SquareDifferenceVect_red);
            [valDiff_green idx_green_1] = min(SquareDifferenceVect_green);
            redSpot(frame, :) = [centroids_red(idx_red_1, 1), centroids_red(idx_red_1, 2)];
            greenSpot(frame, :) = [centroids_green(idx_green_1, 1), centroids_green(idx_green_1, 2)];
            red_green_difference = ((redSpot(frame, 1) - greenSpot(frame, 1))^2) + ((redSpot(frame, 2) - greenSpot(frame, 2))^2)
        
        else
            SquareDifferenceVect_red = (centroids_red(:,1) - centroids_red_last_frame_previous_vid(:,1)).^2 + (centroids_red(:,2) - centroids_red_last_frame_previous_vid(:,2)).^2;
            SquareDifferenceVect_green = (centroids_green(:,1) - centroids_green_last_frame_previous_vid(:,1)).^2 + (centroids_green(:,2) - centroids_green_last_frame_previous_vid(:,2)).^2;
            [valDiff_red, idx_red_1] = min(SquareDifferenceVect_red);
            [valDiff_green, idx_green_1] = min(SquareDifferenceVect_green);
            redSpot(frame, :) = [centroids_red(idx_red_1, 1), centroids_red(idx_red_1, 2)];
            greenSpot(frame, :) = [centroids_green(idx_green_1, 1), centroids_green(idx_green_1, 2)];
        end
        %}
    elseif (frame == first_frame) && ((no_red_led == 1) || (no_green_led ==1))
        first_frame = first_frame + 1;  %******************************
    else
        if (~isempty(centroids_red))
            SquareDifferenceVect_red = (centroids_red(:,1) - redSpot(frame-1, 1)).^2 + (centroids_red(:,2) - redSpot(frame-1, 2)).^2;
            [valDiff_red, idx_red] = min(SquareDifferenceVect_red);
            if (valDiff_red < difference_threshold_red)
                redSpot(frame, :) = [centroids_red(idx_red, 1), centroids_red(idx_red, 2)];
                SquareDifference_red = (redSpot(frame,1) - redSpot(frame-1, 1)).^2 + (redSpot(frame,2) - redSpot(frame-1, 2)).^2;
            else
                redSpot(frame, :) = redSpot(frame-1, :);
                SquareDifference_red = (redSpot(frame,1) - redSpot(frame-1, 1)).^2 + (redSpot(frame,2) - redSpot(frame-1, 2)).^2;
            end
        else
            redSpot(frame, :) = redSpot(frame-1, :);
            SquareDifferenceVect_red = (redSpot(frame,1) - redSpot(frame-1, 1)).^2 + (redSpot(frame,2) - redSpot(frame-1, 2)).^2;
        end
        
        if (~isempty(centroids_green))
            SquareDifferenceVect_green = (centroids_green(:,1) - greenSpot(frame-1, 1)).^2 + (centroids_green(:,2) - greenSpot(frame-1, 2)).^2;
            [valDiff_green, idx_green] = min(SquareDifferenceVect_green);
            greenSpot(frame, :) = [centroids_green(idx_green, 1), centroids_green(idx_green, 2)];
            if (valDiff_green < difference_threshold_green)
                greenSpot(frame, :) = [centroids_green(idx_green, 1), centroids_green(idx_green, 2)];
                SquareDifference_green = (greenSpot(frame,1) - greenSpot(frame-1, 1)).^2 + (greenSpot(frame,2) - greenSpot(frame-1, 2)).^2;
            else
                greenSpot(frame, :) = greenSpot(frame-1, :);
                SquareDifference_green = (greenSpot(frame,1) - greenSpot(frame-1, 1)).^2 + (greenSpot(frame,2) - greenSpot(frame-1, 2)).^2;
            end
        else
            greenSpot(frame, :) = greenSpot(frame-1, :);
            SquareDifference_green = (greenSpot(frame,1) - greenSpot(frame-1, 1)).^2 + (greenSpot(frame,2) - greenSpot(frame-1, 2)).^2;
        end
        
        %if ((SquareDifference_red == 0) && (SquareDifference_green > 5)) || ((SquareDifference_red > 5) && (SquareDifference_green == 0)) || ((SquareDifference_red == 0) && (SquareDifference_green == 0))
        if ((SquareDifference_red == 0) && (SquareDifference_green == 0)) || ((((redSpot(frame, 1) - greenSpot(frame, 1))^2) + ((redSpot(frame, 2) - greenSpot(frame, 2))^2)) > red_green_difference + 120)  %was 120
            
            imshow(thisFrame);
            hold on
            if (~isempty(position_track))
                plot(position_track(:,1), position_track(:,2), 'w')
            end
            hold on
            if (~isempty(centroids_red))
                plot(centroids_red(:,1), centroids_red(:,2), 'r*')
            else
                centroids_red = centroids_red_last_frame_previous_vid;
                plot(centroids_red(:,1), centroids_red(:,2), 'r*')
            end
            hold on
            if (~isempty(centroids_green))
                plot(centroids_green(:,1), centroids_green(:,2), 'g*')
            else
                centroids_green = centroids_green_last_frame_previous_vid;
                plot(centroids_green(:,1), centroids_green(:,2), 'g*')
            end
            display('Choose red spot (use right-click)')
            [x_red, y_red] = getpts     %****************************
            display('Choose green spot (use right-click)')
            [x_green, y_green] = getpts  %***************************
            SquareDifferenceVect_red = [(centroids_red(:,1) - x_red).^2 + (centroids_red(:,2) - y_red).^2];
            SquareDifferenceVect_green = [(centroids_green(:,1) - x_green).^2 + (centroids_green(:,2) - y_green).^2];
            [valDiff_red idx_red_1] = min(SquareDifferenceVect_red);
            [valDiff_green idx_green_1] = min(SquareDifferenceVect_green);
            redSpot(frame, :) = [centroids_red(idx_red_1, 1), centroids_red(idx_red_1, 2)];
            greenSpot(frame, :) = [centroids_green(idx_green_1, 1), centroids_green(idx_green_1, 2)];
        end
    end
    
    if (no_red_led == 0) && (no_green_led == 0)
        
        position_track = [position_track;(redSpot(frame, :) + greenSpot(frame, :))/2];
        
        cla
        imshow(thisFrame);
        imshow(thisFrame, 'Border', 'tight');
        hold on
        plot(position_track(:,1), position_track(:,2), 'y')
        hold on
        plot(redSpot(frame, 1), redSpot(frame, 2), 'r*')
        hold on
        plot(greenSpot(frame, 1), greenSpot(frame, 2), 'g*')
        
        set(gca, 'position', [0 0 1 1], 'units', 'normalized')
        
        %                  filename = [sprintf('%04d',frame) '.tiff'];    %to save video
        %                  saveas(gcf, filename);
        
        %     hold on
        %     line(centroids(:,1),centroids(:,2))
        
        HD_180_deg(frame) = atan2d((greenSpot(frame, 1)-redSpot(frame, 1)), (greenSpot(frame, 2)-redSpot(frame, 2)));
    end
    
end

HD_180_deg_total = [HD_180_ deg_total HD_180_deg];


%{
    num = vid.FramesAcquired;
    frame=getsnapshot(vid);                         % Frame extraction
    framer = frame;
    frameg= framer;
    
    %Red extraction
    diffFrameRed = imsubtract(framer(:,:,1), rgb2gray(framer)); % Get red component of the image
    diffFrameRed = medfilt2(diffFrameRed, [3 3]);   % Filter out the noise by using median filter
    binFrameRed = imbinarize(diffFrameRed, redThresh); % Convert the image into binary image with the red objects as white
    
    %CC_red = bwconncomp(binFrameRed);
    %S_red = regionprops(CC_red, 'Centroid');
    %centroids_red = cat(1, S_red.Centroid);
    
    %Green extraction
    diffFrameGreen = imsubtract(frameg(:,:,2), rgb2gray(frameg)); % Get green component of the image
    diffFrameGreen = medfilt2(diffFrameGreen, [3 3]); % Filter out the noise by using median filter
    binFrameGreen = imbinarize(diffFrameGreen, greenThresh); % Convert the image into binary image with the green objects as white
    
    %CC_green = bwconncomp(binFrameGreen);
    %S_green = regionprops(CC_green, 'Centroid');
    %centroids_green = cat(1, S_green.Centroid);
    
    %live tracking/Identifying objects
    %RED
    diffFrameRed = bwareaopen(binFrameRed,75);    %eliminates any object smaller than the specified pxls
    rlabel = bwlabel(diffFrameRed,8);               %Labels connectected components
    objR = regionprops(rlabel,'BoundingBox','Centroid');
    imshow(frame)
    hold on
    
    %Green
    diffFrameGreen = bwareaopen(binFrameGreen,50);    %eliminates any object smaller than the specified pxls
    glabel = bwlabel(diffFrameGreen,8);               %Labels connectected components
    objG = regionprops(glabel,'BoundingBox','Centroid');
    imshow(frame)
    hold on
    
    if num ~=0
        %Red object creation around the centroid and data collection if red present
        if ~isempty(objR)
            for objectR = 1:length(objR)
                bcr = objR(objectR).BoundingBox;
                cr = objR(objectR).Centroid;
                %Boundary condition enforcement
                if((objR(objectR).Centroid(1)<= 10) || (objR(objectR).Centroid(2)<= 10) || (objR(objectR).Centroid(1)>= X_resolution-10) || (objR(objectR).Centroid(2)>= Y_resolution-10))
                    REDco(num, 2*objectR-1) = REDco(num-1, 2*objectR-1);
                    REDco(num, 2*objectR) = REDco(num-1, 2*objectR);
                else
                    REDco(num, 2*objectR-1) = objR(objectR).Centroid(1);
                    REDco(num, 2*objectR) = objR(objectR).Centroid(2);
                    rectangle('Position',bcr,'EdgeColor','r','LineWidth',2)
                    plot(cr(1),cr(2),'-m+')
                end
                if abs(REDco(num, 2*objectR-1)-(X_resolution/2)<= CentCo(1,1)) && abs(REDco(num, 2*objectR)-(Y_resolution) <= CentCo(1,2))
                    CentCo(1,1) = REDco(num, 2*objectR-1);
                    CentCo(1,2) = REDco(num, 2*objectR);
                end
            end
            %track the point nearest to the center (to avoid wall reflection)
            %if REDco(num, 1) <
        end
        %Green Object creation around the centroid and data collection if green is present
        if ~isempty(objG)
            for objectG = 1:length(objG)
                bcg = objG(objectG).BoundingBox;
                cg = objG(objectG).Centroid;
                if((objG(objectG).Centroid(1)<= 10) || (objG(objectG).Centroid(2)<= 10) || (objG(objectG).Centroid(1)>= X_resolution-10) || (objG(objectG).Centroid(2)>= Y_resolution-10))
                    GREENco(num, 2*objectG-1) = GREENco(num-1, 2*objectG-1);
                    GREENco(num, 2*objectG) = GREENco(num-1, 2*objectG);
                else
                    GREENco(num, 2*objectG-1) = objG(objectG).Centroid(1);
                    GREENco(num, 2*objectG) = objG(objectG).Centroid(2);
                end
                rectangle('Position',bcg,'EdgeColor','g','LineWidth',2)
                plot(cg(1),cg(2),'-m+')
            end
        end
    end
    
%     HD_180_deg(frame) = atan2d((greenSpot(frame, 1)-redSpot(frame, 1)), (greenSpot(frame, 2)-redSpot(frame, 2)));
    plot(GREENco(:,1),GREENco(:,2),'g',REDco(:,1),REDco(:,2),'r')
    %drawnow
    hold off
%}

stop(vid);
flushdata(vid);

figure(2)
plot(REDco(:,1),-REDco(:,2),'r',GREENco(:,1),-GREENco(:,2),'g')
