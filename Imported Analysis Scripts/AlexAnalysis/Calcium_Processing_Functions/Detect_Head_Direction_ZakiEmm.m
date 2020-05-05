% Number of frames to inspect
FrameNum = 100;

FrameNum = 100;  % Number of frames to inspect                                    
REDco = NaN(FrameNum, 5); %Red Centroid coordinate storage
GREENco = NaN(FrameNum,5);%Green Centroid Coordinate storage
HEADco = zeros(FrameNum,5); %Head postion coordinate storage
HDdeg = zeros(FrameNum, 1);
redThresh = 0.007;   % 0.01 * 0.0039
greenThresh = 0.007;  % 0.01 //best: 0.007
difference_threshold_red = 0.9;  % 250   //best: 350
difference_threshold_green = 0.4;  % 260  //best: 400  700
ExtentLim = 0.4;                %Extent limit : ratio of bounding box area to full area
EccLim = 0.9;                   %Eccentricity limit: roundness factor, 1=line, 0=circle 
GreenIndex = 0;
vid = [];
FrameInterval = 1;
frame = 1;

for ind_vid = 1:9
    % Open the VT1.avi demo movie that ships with MATLAB.
    % First get the folder that it lives in.
    folder = ['X:\Emmanuel\Data\HDtracking'];
    movieFullFileName = ['behavCam1.avi'];
    folder = fileparts(which(['behavCam',num2str(ind_vid),'.avi'])); % Determine where video folder is.
    movieFullFileName = fullfile(folder, ['behavCam',num2str(ind_vid),'.avi']);
    % Check to see that it exists.
    if ~exist(movieFullFileName, 'file')
        strErrorMessage = sprintf('File not found:\n%s\nYou can choose a new one, or cancel', movieFullFileName);
        response = questdlg(strErrorMessage, 'File not found', 'OK - choose a new movie.', 'Cancel', 'OK - choose a new movie.');
        if strcmpi(response, 'OK - choose a new movie.')
            [baseFileName, folderName, FilterIndex] = uigetfile('*.avi');
            if ~isequal(baseFileName, 0)
                movieFullFileName = fullfile(folderName, baseFileName);
            else
                return;
            end
        else
            return;
        end
    end
    
    vid = VideoReader(movieFullFileName);
    % Determine how many frames there are.
    FrameNum = vid.NumberOfFrames;
    X_resolution = vid.Height;
    Y_resolution = vid.Width;
    frame =1;
    
    while(frame <= FrameNum)               %# of frames being processed
        
        num = frame;
        %frame=getsnapshot(vid);                         % Frame extraction
        framer = read(vid, frame);
        frameg= read(vid,frame);
        
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
        diffFrameRed = bwareaopen(binFrameRed,40);    %eliminates any object smaller than the specified pxls
        CC_red = bwconncomp(diffFrameRed);
        objR = regionprops(diffFrameRed,'BoundingBox','Centroid', 'Eccentricity', 'Extent'); %Labels connectected components with boxed and label/identify centroid
        idr = find([objR.Eccentricity] < EccLim & [objR.Extent] > ExtentLim);                % Identifies all objects that do not fit the required description
        framer = ismember(labelmatrix(CC_red), idr);                                         % Eliminate all objects that do not comply to the previous line's conditions
        targetR = regionprops(framer,'BoundingBox', 'Centroid');                             % Redefine image with only potentially wanted objects
        imshow(frame)
        hold on
        
        %Green
        diffFrameGreen = bwareaopen(binFrameGreen,40);    %eliminates any object smaller than the specified pxls
        CC_Green = bwconncomp(diffFrameGreen);
        objG = regionprops(diffFrameGreen,'BoundingBox','Centroid','Eccentricity', 'Extent');%Labels connectected components with boxed and label/identify centroid
        idg = find([objG.Eccentricity] < EccLim & [objG.Extent] > ExtentLim);   % Identifies all objects that do not fit the required description
        frameg = ismember(labelmatrix(CC_Green), idg);                          % Eliminate all objects that do not comply to the previous line's conditions
        targetG = regionprops(frameg,'BoundingBox', 'Centroid');                % Redefine image with only potentially wanted objects
        imshow(frame)
        hold on
        
        if num ~=0 && ~isempty(targetR)&& ~isempty(targetG)
            %Red object creation around the centroid and data collection if red present
            for objectR = 1:length(targetR)    %for every red object
                for objectG = 1:length(targetG)
                    %find any two points that are within a respectable distance from each other in both the x and y direction
                    if(abs(targetR(objectR).Centroid(1)-targetG(objectG).Centroid(1))<40 && abs(targetR(objectR).Centroid(2)-targetG(objectG).Centroid(2))<30)
                        GreenIndex = objectG;                                      %Green point indices
                        bcr = targetR(objectR).BoundingBox;
                        cr = targetR(objectR).Centroid;
                        REDco(num, 1) = targetR(objectR).Centroid(1);               %Centroid X position coordinate
                        REDco(num, 2) = targetR(objectR).Centroid(2);               %Centroid Y position Coordinate
                        rectangle('Position',bcr,'EdgeColor','r','LineWidth',2)
                        plot(cr(1),cr(2),'-m+')
                    end
                end
            end
            %Green Object creation around the centroid and data collection if green is present
            if GreenIndex ~=0
                bcg = targetG(GreenIndex).BoundingBox;
                cg = targetG(GreenIndex).Centroid;
                GREENco(num, 2*GreenIndex-1) = targetG(GreenIndex).Centroid(1); %Centroid X position coordinate
                GREENco(num, 2*GreenIndex) = targetG(GreenIndex).Centroid(2);   %Centroid Y position Coordinate
                rectangle('Position',bcg,'EdgeColor','g','LineWidth',2)
                plot(cg(1),cg(2),'-m+')
                GreenIndex = 0;
            end
        end
        HEADco(:,1) = abs(GREENco(:,1)+ REDco(:,1))/2;  %Head X position Coordinates
        HEADco(:,2) = abs(GREENco(:,2) + REDco(:,2))/2; %Head Y position Coordinates
        %HDdeg(:,1) = atan2d(GREENco(:,2)-REDco(:,2) , GREENco(:,1) - REDco(:,1));  %Head orientation angle LIVE
        plot(GREENco(:,1),GREENco(:,2),'g.',REDco(:,1),REDco(:,2),'r.',HEADco(:,1),HEADco(:,2),'w')    %Live plot
        hold off
        frame = frame + int16(FrameNum/FrameInterval);
    end
end