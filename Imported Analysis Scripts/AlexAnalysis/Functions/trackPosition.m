function [pos hd] = trackPosition(doPath)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Head Tracking Script    %
    % Author: Emmanuel Wilson %
    % This Script is designed to locate and track green and red LED's found
    % on a mouses' head and extrapolate the head location and orientation for
    % post processing. It will output an LED tracking plot/head position on top of the
    % behavioural video, a plot of overall tracking and a histogram of
    % non-tracked/missed frames. Make sure that all your behavioral ".avi" videos 
    % are in the current path and run!
    %Please download: https://www.mathworks.com/matlabcentral/fileexchange/1039-hline-and-vline
    %and place files in your file path for horizontal/vertical line segments in your plots. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % clc
    % close all
    % clear

    folder = dir(doPath);     %Looks at current file location

    totalFrames = 0;        %Total number of frames observed (including previous video itterations)
    TimestampMin = 10000000;%Arbitrary VERY large number for min value storage, must be larger than datenum stamp in folder
    Timestamp = 0;          %TimeStamp record for indexing
    o = NaN(1,length(folder)); %index movie order

    % find avi and dat files
    aviFiles = dir([doPath '\*.avi']);
    filePrefix = 'behavCam';

    %Find total number of frames and organize vidoes by their timestamps 
    for ind_vid = 1: length(folder)                         %Going through every object in the folder
        if (~isempty(strfind(folder(ind_vid).name,filePrefix)))
            movieFullFileName =  folder(ind_vid).name;      %Movie file name
            video = VideoReader([doPath '/' movieFullFileName]);   %Read video file              
            totalFrames = totalFrames + video.NumberOfFrames;    %Determine how many frames in file
            for t = 1 : length(folder)
                if (~isempty(strfind(folder(t).name,filePrefix)) && folder(t).datenum < TimestampMin && folder(t).datenum > Timestamp)
                    TimestampMin = folder(t).datenum;       %Save as minimum value
                    o(1,ind_vid) = t;                       %Save Index #
                end
            end
            Timestamp = TimestampMin;                       %Set minimum time boundary
            TimestampMin = 10*TimestampMin;                 %Reset minimum value
        end               
    end

    videoObject = [];
    HEADco = zeros(totalFrames,5);     %Head coordinates storage. Each instance is a frame
    PLOTdata = NaN(totalFrames, 5);    %Modified head coordinates storage
    REDco = NaN(totalFrames, 5);       %Red Centroid coordinate storage Each instance is a frame
    GREENco = NaN(totalFrames,5);      %Green Centroid Coordinate storage Each instance is a frame
    HDdeg = zeros(totalFrames, 1);     %Head direction/orientation storage
    SINKdata = HEADco;                 %Head coordinates used for sinking with the miniscope
    tempHDdeg = zeros(totalFrames, 1); % temporary Head direction/orientation storage
    redThresh = 0.75;          % Must be a value between 1 and 0. The lower the number the more sensitve it is
    greenThresh = 0.75;        % Must be a value between 1 and 0. The lower the number the more sensitve it is
    ExtentLim = 0.2;            %Extent limit : ratio of bounding box area to full area
    EccLim = 0.9;               %Eccentricity limit: roundness factor, 1=line, 0=circle
    GreenIndex = 0;             %Flag indicating if green object (LED) has been located
    LEDxDist = 20;              %Max distance between LED's in the x axis in pixels
    LEDyDist = 15;              %Max distance between LED's in the y axis in pixels
    frame  = 1;                 %Frame # of current video
    CurrentFrame = frame;

    for it = 1 : length(o)                 %Going through every object in the folder
        ind_vid = o(1,it);
        if ~isnan(ind_vid) && (folder(ind_vid).isdir == 0)                     %ignore if file is not a video
            movieFullFileName =  folder(ind_vid).name;      %Movie file name
            videoObject = VideoReader([doPath '/' movieFullFileName]);   %Read video file      
            numberOfFrames = videoObject.NumberOfFrames;    %Determine how many frames in current video file   
            
            
            v = VideoReader([doPath '/' movieFullFileName]);             %Read video file a gain ( hasFrame() function won't work otherwise)
           
            % Loop through the movie.
%             meanFrame = read(videoObject, 1);
%             frame = 2;
%             while hasFrame(v) && frame<= numberOfFrames
%                 meanFrame = nansum(cat(4,read(videoObject, frame),meanFrame),4);
%                 frame = frame + 1;
%             end
%             
%             meanFrame = meanFrame ./ (frame-1);
            
            fprintf(['\n\t\tTracking Frames (Clip ' num2str(it-2) '):  '])
            doBack = 0;
            tic
            while hasFrame(v) && frame<= numberOfFrames
%                 ns = [num2str(round(1000.*frame./numberOfFrames)./10) ' %%'];
%                 fprintf([repmat('\b',[1 doBack]) ns]);
%                 doBack = length(ns)-1;
                
                
                thisFrame = read(videoObject, frame);       % Extract the frame from the movie structure.
                framer = thisFrame;                         % frame for red analysis
                frameg= framer;                             % frame for green analysis
                %Red extraction
                diffFrameRed = imsubtract(framer(:,:,1), [framer(:,:,2)+framer(:,:,3)]./2); % Get red component of the image
                diffFrameRed = medfilt2(diffFrameRed, [3 3]);   % Filter out the noise by using median filter
                binFrameRed = (diffFrameRed > 10);
                
                %Green extraction
                diffFrameGreen = imsubtract(framer(:,:,2), [framer(:,:,1)+framer(:,:,3)]./2); % Get green component of the image
                diffFrameGreen = medfilt2(diffFrameGreen, [3 3]); % Filter out the noise by using median filter
                binFrameGreen = (diffFrameGreen > 3);
                
                %live tracking/Identifying objects
                %RED
                diffFrameRed = bwareaopen(binFrameRed,10);    %eliminates any object smaller than the specified pxls
                CC_red = bwconncomp(diffFrameRed);
                objR = regionprops(diffFrameRed,'BoundingBox','Centroid', 'Eccentricity', 'Extent'); %Labels connectected components with boxed and label/identify centroid
                idr = find([objR.Eccentricity] < EccLim & [objR.Extent] > ExtentLim);                % Identifies all objects that do not fit the required description
                framer = ismember(labelmatrix(CC_red), idr);                                         % Eliminate all objects that do not comply to the previous line's conditions
                targetR = regionprops(framer,'BoundingBox', 'Centroid');                             % Redefine image with only potentially wanted objects
%                 imshow(thisFrame)
%                 hold on

                %Green
                diffFrameGreen = bwareaopen(binFrameGreen,2);    %eliminates any object smaller than the specified pxls
                CC_Green = bwconncomp(diffFrameGreen);
                objG = regionprops(diffFrameGreen,'BoundingBox','Centroid','Eccentricity', 'Extent');%Labels connectected components with boxed and label/identify centroid
                idg = find([objG.Eccentricity] < EccLim & [objG.Extent] > ExtentLim);   % Identifies all objects that do not fit the required description
                frameg = ismember(labelmatrix(CC_Green), idg);                          % Eliminate all objects that do not comply to the previous line's conditions
                targetG = regionprops(frameg,'BoundingBox', 'Centroid');                % Redefine image with only potentially wanted objects
%                 imshow(thisFrame)
%                 hold on

                if frame ~=0 && ~isempty(targetR)&& ~isempty(targetG)                   %skip frame unless green and red can be found
                    %Red object creation around the centroid and data collection if red present
                    for objectR = 1:length(targetR)                                     %for every red object...
                        for objectG = 1:length(targetG)                                 %for every green object...
                            %find any two points that are within a respectable distance from each other in both the x and y direction
                            if(abs(targetR(objectR).Centroid(1)-targetG(objectG).Centroid(1))<LEDxDist && abs(targetR(objectR).Centroid(2)-targetG(objectG).Centroid(2))<LEDyDist)
                                GreenIndex = objectG;                                    %Green point indices
                                %Centroid info extraction
                                REDco(CurrentFrame, 1) = targetR(objectR).Centroid(1);    %Centroid X position coordinate
                                REDco(CurrentFrame, 2) = targetR(objectR).Centroid(2);    %Centroid Y position Coordinate
%                                 rectangle('Position',bcr,'EdgeColor','r','LineWidth',2)  %creating visual representation of bounding box
%                                 plot(cr(1),cr(2),'-m+')                                  %plot on current frame
                            end
                        end
                    end
                    %Green Object creation around the centroid and data collection if green is present
                    if GreenIndex ~=0
                           %Centroid info extraction
                        GREENco(CurrentFrame, 2*GreenIndex-1) = targetG(GreenIndex).Centroid(1); %Centroid X position coordinate
                        GREENco(CurrentFrame, 2*GreenIndex) = targetG(GreenIndex).Centroid(2);   %Centroid Y position Coordinate
%                         rectangle('Position',bcg,'EdgeColor','g','LineWidth',2)         %creating visual representation of bounding box
%                         plot(cg(1),cg(2),'-m+')                                         %plot on current frame
                        GreenIndex = 0;                                                 %reset green index
                    end
                end
                
                HEADco(:,1) = abs(GREENco(:,1)+ REDco(:,1))/2;  %Head X position Coordinates
                HEADco(:,2) = abs(GREENco(:,2) + REDco(:,2))/2; %Head Y position Coordinates
                tempHDdeg(:,1) = 180 + atan2d(GREENco(:,2)-REDco(:,2) , GREENco(:,1) - REDco(:,1));  %Head orientation angle
%                 plot(GREENco(:,1),GREENco(:,2),'g.',REDco(:,1),REDco(:,2),'r.',HEADco(:,1),HEADco(:,2),'w')    %Live plot
%                 hold off
                frame = frame + 1;
                CurrentFrame = CurrentFrame + 1;
            end
            toc
            fprintf('\b')
        end
        frame = 1;
    end
    
    smoothedPos = interpNaNs(HEADco(:,1:2));
    isJump = sqrt(sum(diff(smoothedPos(:,1:2)).^2,2))>6;
    for i = -7:1:7
        smoothedPos(circshift(isJump,[i 0]),:) = nan;
    end
    pos = interpNaNs(smoothedPos(:,1:2))';
    hd = tempHDdeg(:,1)';
    
    
    close all
    figure(1)
    set(gcf,'position',[50 50 300 300])
    plot(pos(:,1),pos(:,2))
    axis equal
    drawnow
    inds = find(ismember(doPath,'/'));
    outP = doPath;
    outP(inds) = '_';
    outP = outP(inds+1:end);
    saveFig(gcf,['Plots/PositionTracking/' outP],'tiff');
end



































