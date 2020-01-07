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

folder = dir(pwd);     %Looks at current file location

totalFrames = 0;        %Total number of frames observed (including previous video itterations)
TimestampMin = 10000000;%Arbitrary VERY large number for min value storage, must be larger than datenum stamp in folder
Timestamp = 0;          %TimeStamp record for indexing
o = NaN(1,length(folder)); %index movie order

% find avi and dat files
aviFiles = dir([folder.name '\*.avi']);
filePrefix = 'behavCam';
count = 0;

%Sorting through file timestamps and marking the index's for future reference
for i = 1 : ms.numFiles
    for j = 1 : length(folder)
        if(folder(j).datenum == Timestamp)
            count = count +1;
            anomilynames(1,count) = j;
        end
        if (~isempty(strfind(folder(j).name,filePrefix)) && folder(j).datenum < Timestampmin && folder(j).datenum > Timestamp)
            o(1,i) = j;                         %Store file location within folder
            Timestampmin = folder(j).datenum;   %Set minimum time to lowest value found thought the loop
            %                 prevname = folder(j).name;
        end
        if count > 2
            if folder(anomilynames(1,2)).bytes < folder(anomilynames(1,3)).bytes
                o(1,i) = anomilynames(1,3);                         %Store file location within folder
                Timestampmin = folder(j).datenum;   %Set minimum time to lowest value found thought the loop
            end
        end
    end
    Timestamp = Timestampmin;                   %Reset loop boundaries
    Timestampmin = Timestampmin*100;             %Reset loop boundaries
    count = 0;
    
end
if length(o)>1 && folder(o(1,length(o(1,:)))).bytes > folder(o(1,length(o(1,:))-1)).bytes
    temp = o(1,length(o(1,:)));
    o(1,length(o(1,:)))= o(1,length(o(1,:))-1);
    o(1,length(o(1,:))-1)= temp;
end

videoObject = [];
HEADco = zeros(totalFrames,5);     %Head coordinates storage. Each instance is a frame
PLOTdata = NaN(totalFrames, 5);    %Modified head coordinates storage
REDco = NaN(totalFrames, 5);       %Red Centroid coordinate storage Each instance is a frame
GREENco = NaN(totalFrames,5);      %Green Centroid Coordinate storage Each instance is a frame
HDdeg = zeros(totalFrames, 1);     %Head direction/orientation storage
SINKdata = HEADco;                 %Head coordinates used for sinking with the miniscope
tempHDdeg = zeros(totalFrames, 1); % temporary Head direction/orientation storage
redThresh = 0.007;          % Must be a value between 1 and 0. The lower the number the more sensitve it is
greenThresh = 0.007;        % Must be a value between 1 and 0. The lower the number the more sensitve it is
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
        videoObject = VideoReader(movieFullFileName);   %Read video file
        numberOfFrames = videoObject.NumberOfFrames;    %Determine how many frames in current video file
        v = VideoReader(movieFullFileName);             %Read video file a gain ( hasFrame() function won't work otherwise)
        % Loop through the movie.
        while hasFrame(v) && frame<= numberOfFrames
            thisFrame = read(videoObject, frame);       % Extract the frame from the movie structure.
            framer = thisFrame;                         % frame for red analysis
            frameg= framer;                             % frame for green analysis
            %Red extraction
            diffFrameRed = imsubtract(framer(:,:,1), rgb2gray(framer)); % Get red component of the image
            diffFrameRed = medfilt2(diffFrameRed, [3 3]);   % Filter out the noise by using median filter
            binFrameRed = imbinarize(diffFrameRed, redThresh); % Convert the image into binary image with the red objects as white
            
            %Green extraction
            diffFrameGreen = imsubtract(frameg(:,:,2), rgb2gray(frameg)); % Get green component of the image
            diffFrameGreen = medfilt2(diffFrameGreen, [3 3]); % Filter out the noise by using median filter
            binFrameGreen = imbinarize(diffFrameGreen, greenThresh); % Convert the image into binary image with the green objects as white
            
            %live tracking/Identifying objects
            %RED
            diffFrameRed = bwareaopen(binFrameRed,10);    %eliminates any object smaller than the specified pxls
            CC_red = bwconncomp(diffFrameRed);
            objR = regionprops(diffFrameRed,'BoundingBox','Centroid', 'Eccentricity', 'Extent'); %Labels connectected components with boxed and label/identify centroid
            idr = find([objR.Eccentricity] < EccLim & [objR.Extent] > ExtentLim);                % Identifies all objects that do not fit the required description
            framer = ismember(labelmatrix(CC_red), idr);                                         % Eliminate all objects that do not comply to the previous line's conditions
            targetR = regionprops(framer,'BoundingBox', 'Centroid');                             % Redefine image with only potentially wanted objects
            imshow(thisFrame)
            hold on
            
            %Green
            diffFrameGreen = bwareaopen(binFrameGreen,2);    %eliminates any object smaller than the specified pxls
            CC_Green = bwconncomp(diffFrameGreen);
            objG = regionprops(diffFrameGreen,'BoundingBox','Centroid','Eccentricity', 'Extent');%Labels connectected components with boxed and label/identify centroid
            idg = find([objG.Eccentricity] < EccLim & [objG.Extent] > ExtentLim);   % Identifies all objects that do not fit the required description
            frameg = ismember(labelmatrix(CC_Green), idg);                          % Eliminate all objects that do not comply to the previous line's conditions
            targetG = regionprops(frameg,'BoundingBox', 'Centroid');                % Redefine image with only potentially wanted objects
            imshow(thisFrame)
            hold on
            
            if frame ~=0 && ~isempty(targetR)&& ~isempty(targetG)                   %skip frame unless green and red can be found
                %Red object creation around the centroid and data collection if red present
                for objectR = 1:length(targetR)                                     %for every red object...
                    for objectG = 1:length(targetG)                                 %for every green object...
                        %find any two points that are within a respectable distance from each other in both the x and y direction
                        if(abs(targetR(objectR).Centroid(1)-targetG(objectG).Centroid(1))<LEDxDist && abs(targetR(objectR).Centroid(2)-targetG(objectG).Centroid(2))<LEDyDist)
                            GreenIndex = objectG;                                    %Green point indices
                            bcr = targetR(objectR).BoundingBox;                      %bounding box info extraction
                            cr = targetR(objectR).Centroid;                          %Centroid info extraction
                            REDco(CurrentFrame, 1) = targetR(objectR).Centroid(1);    %Centroid X position coordinate
                            REDco(CurrentFrame, 2) = targetR(objectR).Centroid(2);    %Centroid Y position Coordinate
                            rectangle('Position',bcr,'EdgeColor','r','LineWidth',2)  %creating visual representation of bounding box
                            plot(cr(1),cr(2),'-m+')                                  %plot on current frame
                        end
                    end
                end
                %Green Object creation around the centroid and data collection if green is present
                if GreenIndex ~=0
                    bcg = targetG(GreenIndex).BoundingBox;                          %bounding box info extraction
                    cg = targetG(GreenIndex).Centroid;                              %Centroid info extraction
                    GREENco(CurrentFrame, 2*GreenIndex-1) = targetG(GreenIndex).Centroid(1); %Centroid X position coordinate
                    GREENco(CurrentFrame, 2*GreenIndex) = targetG(GreenIndex).Centroid(2);   %Centroid Y position Coordinate
                    rectangle('Position',bcg,'EdgeColor','g','LineWidth',2)         %creating visual representation of bounding box
                    plot(cg(1),cg(2),'-m+')                                         %plot on current frame
                    GreenIndex = 0;                                                 %reset green index
                end
            end
            HEADco(:,1) = abs(GREENco(:,1)+ REDco(:,1))/2;  %Head X position Coordinates
            HEADco(:,2) = abs(GREENco(:,2) + REDco(:,2))/2; %Head Y position Coordinates
            tempHDdeg(:,1) = 180 + atan2d(GREENco(:,2)-REDco(:,2) , GREENco(:,1) - REDco(:,1));  %Head orientation angle
            plot(GREENco(:,1),GREENco(:,2),'g.',REDco(:,1),REDco(:,2),'r.',HEADco(:,1),HEADco(:,2),'w')    %Live plot
            hold off
            frame = frame + 1;
            CurrentFrame = CurrentFrame + 1;
        end
    end
    frame = 1;
end

%% Post-post Processing, interpolation and error stats
%Only interpolate 5 missed frames in a row or lower
y = 1;      %good frame index
nany = 0;   %#of bad frames (NaN) in a row
Er = 0;     %Error/bad frame flag
ErHist = NaN(totalFrames,1);        %Histogram data storage
nancount = 0;                       %total Error/NaN count
nangcount = 0;                      %Error/NaN count < 5
for i = 1 : length(HEADco)                                          %For the whole data set
    if (~isnan(HEADco(i,1)) && Er == 0) || (~isnan(HEADco(i,1)) && (Er == 1 && nany > 5))   %if HEADco is not NaN and not an error or if there were more than 5 errors in a row then ...
        PLOTdata(y,:) = HEADco(i,:);                                %transfer head postion to plot data
        SINKdata(i,:) = HEADco(i,:);                                %transfer head postion to sink data
        HDdeg(i,:) = tempHDdeg(i,:);                                %transfer orientation to HDdeg
        y = y + 1;
        if Er == 1
            ErHist(i,1) = nany;
        end
        nany = 0;                                                   %reset flag
        Er = 0;                                                     %reset flag
    elseif ~isnan(HEADco(i,1)) && Er == 1 && nany <=5               %Interpolate error/empty frames that are 5 frames in a row or less
        intrpX = (HEADco(i,1) - HEADco(i-(nany+1),1))/(nany+1);     %Create even spacing between interpolated points in the x direction
        intrpY = (HEADco(i,2)- HEADco(i-(nany+1),2))/(nany+1);      %Create even spacing between interpolated points in the x direction
        degIn =  (tempHDdeg(i,1) - tempHDdeg(i-(nany+1),1))/(nany+1);%Create even spacing between interpolated points for orientation
        PLOTdata(y,1) = HEADco(i-(nany+1),1)+ intrpX;               %Interpolate incriment in X
        PLOTdata(y,2) = HEADco(i-(nany+1),2)+ intrpY;               %Interpolate incriment in Y
        SINKdata(i-nany,1) = PLOTdata(y,1);                         %Interpolated data in SINK
        SINKdata(i-nany,2) = PLOTdata(y,2);                         %Interpolated data in SINK
        HDdeg(i,1) = tempHDdeg(i-(nany+1),1)+degIn;                 %Interpolate incriment in orientation
        y = y+1;
        if nany >= 2
            for s = 1 : nany-1
                PLOTdata(y,1) = PLOTdata(y-1,1) + intrpX;          %Interpolate incriment in X
                PLOTdata(y,2) = PLOTdata(y-1,2) + intrpY;          %Interpolate incriment in Y
                SINKdata(i-(nany-s),1) = PLOTdata(y,1);            %Interpolated data in SINK
                SINKdata(i-(nany-s),2) = PLOTdata(y,2);            %Interpolated data in SINK
                HDdeg(i,1) = HDdeg(y-1,1)+degIn;                   %Interpolate incriment in orientation
                y = y + 1;
            end
        end
        PLOTdata(y,:) = HEADco(i,:);                    %transfer head postion to plot data
        SINKdata(i,:) = HEADco(i,:);                    %transfer head position to SINK data
        HDdeg(i,:) = tempHDdeg(i,:);                    %transfer orientation to HDdeg
        y = y + 1;
        ErHist(i,1) = nany;                             %Store consecutive NaN values in histogram
        nangcount = nangcount + nany;                   %Increase Good/interpolated NaN frame count
        nany = 0;                                       %reset flag
        Er = 0;                                         %reset flag
    else                                                %if it is an error frame
        nany = nany + 1;                                %incriment # of error/NaN frames in a row
        Er = 1 ;                                        %error flag
        nancount = nancount + 1;                        %Increase NaN count
        SINKdata(i,:) = HEADco(i,:);                    %transfer NaN value from HEADco to SINK
        HDdeg(i,:) = tempHDdeg(i,:);                    %transfer orientation to HDdeg
    end
end

%spatial conversion approximation: cm/pixel
pixConX = 75/length(frameg(1,:));
pixConY = 75/length(frameg(:,1));

%Error stats
ErTOT = ((y-nangcount)/totalFrames)*100;        % Percent of good files pre-interpolation
Ertopf = (y/totalFrames)*100;                   % Percent of good files post-interpolation
%Plot head postion of entire session
figure(2)
plot(PLOTdata(:,1),-PLOTdata(:,2),'b')
title('Trajectory plot')
%Plot error Histogram with error stats
h = figure(3)
histogram(ErHist)
xlabel('Number of consecutive error/NaN frames')
ylabel('Number of instances')
vline(5.5)
stats = char (['Usable frames before interpolation: ', num2str(ErTOT), '%'], ['Usable frames after interpolation: ', num2str(Ertopf), '%']);
text(7,100,stats)
title('Error Report; Frames Missed')

save('HeadTrackingData','SINKdata','HDdeg','pixConX','pixConY')
savefig(h,'ErrorReportHist')
%% Save video
% folder = fileparts(which(strcat('behavCam', num2str(11), '.avi')));  % msCam8   //// extraction videos 1,2 and 3 are msCam 3, 4 and 5 (3 dec)
% workingDir = folder;
%     imageNames = dir(fullfile(workingDir,['images_',num2str(11)],'*.tiff'));
% imageNames = {imageNames.name}';
%
% outputVideo = VideoWriter(fullfile(workingDir, strcat('extraction_test_', num2str(11), '.avi')));
% open(outputVideo)
%
% for ii = 1:length(imageNames)
%    img = imread(fullfile(workingDir,['images_',num2str(11)],imageNames{ii}));
%    writeVideo(outputVideo,img)
% end
%
% close(outputVideo)