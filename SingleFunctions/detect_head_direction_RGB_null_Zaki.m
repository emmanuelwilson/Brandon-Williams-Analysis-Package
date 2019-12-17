% clc
% close all
% clear all
HD_180_deg_total = [];
position_track = [];
no_red_led = 1;
no_green_led = 1;
centroids_red_last_frame_previous_vid = 0;
centroids_green_last_frame_previous_vid = 0;

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
    videoObject = [];
    videoObject = VideoReader(movieFullFileName)
    % Determine how many frames there are.
    numberOfFrames = videoObject.NumberOfFrames;
    vidHeight = videoObject.Height;
    vidWidth = videoObject.Width;
    
    numberOfFramesWritten = 0;
    
    HD = 0;
    HD_180_deg = 0;
    big_spot = 0;
    small_spot = 0;
    error_frame = [];
    num_error_frames = 0;
    
    redThresh = 0.0035;   % 0.01
    greenThresh = 0.001;  % 0.01 //best: 0.007
    difference_threshold_red = 500;  % 250   //best: 350
    difference_threshold_green = 700;  % 260  //best: 400  700
    
    first_frame = 1;
    frame  = first_frame;
    
    % Loop through the movie.
    while (frame <= numberOfFrames)
        % Extract the frame from the movie structure.
        thisFrame = read(videoObject, frame);
        
        diffFrameRed = imsubtract(thisFrame(:,:,1), rgb2gray(thisFrame)); % Get red component of the image
        diffFrameRed = medfilt2(diffFrameRed, [3 3]); % Filter out the noise by using median filter
        binFrameRed = im2bw(diffFrameRed, redThresh); % Convert the image into binary image with the red objects as white
        
        CC_red = bwconncomp(binFrameRed);
        S_red = regionprops(CC_red, 'Centroid');
        
        centroids_red = cat(1, S_red.Centroid);
        
        if ~isempty(centroids_red)
            centroids_red_last_frame_previous_vid = centroids_red;
            if ind_vid == 1
            no_red_led = 0;
            end
        end
        
        %     redPixelIdxList_lengthVect = [];
        %     for i = 1 : length(CC_red.PixelIdxList)
        %         redPixelIdxList_lengthVect(i) = length(CC_red.PixelIdxList{i});
        %     end
        %
        %     [val_red idx_red] = max(redPixelIdxList_lengthVect);
        %     redSpot(frame, :) = [centroids_red(idx_red, 1), centroids_red(idx_red, 2)];
        
        diffFrameGreen = imsubtract(thisFrame(:,:,2), rgb2gray(thisFrame)); % Get green component of the image
        diffFrameGreen = medfilt2(diffFrameGreen, [3 3]); % Filter out the noise by using median filter
        binFrameGreen = im2bw(diffFrameGreen, greenThresh); % Convert the image into binary image with the green objects as white
        
        CC_green = bwconncomp(binFrameGreen);
        S_green = regionprops(CC_green, 'Centroid');
        
        centroids_green = cat(1, S_green.Centroid);
        
        if ~isempty(centroids_green)
            centroids_green_last_frame_previous_vid = centroids_green;
            if ind_vid == 1;
            no_green_led = 0;
            end
        end
        
        %     greenPixelIdxList_lengthVect = [];
        %     for ii = 1 : length(CC_green.PixelIdxList)
        %         greenPixelIdxList_lengthVect(ii) = length(CC_green.PixelIdxList{ii});
        %     end
        %
        %     [val_green idx_green] = max(greenPixelIdxList_lengthVect);
        %     greenSpot(frame, :) = [centroids_green(idx_green, 1), centroids_green(idx_green, 2)];
        
        
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
                SquareDifferenceVect_red = [(centroids_red(:,1) - centroids_red_last_frame_previous_vid(:,1)).^2 + (centroids_red(:,2) - centroids_red_last_frame_previous_vid(:,2)).^2];
                SquareDifferenceVect_green = [(centroids_green(:,1) - centroids_green_last_frame_previous_vid(:,1)).^2 + (centroids_green(:,2) - centroids_green_last_frame_previous_vid(:,2)).^2];
                [valDiff_red idx_red_1] = min(SquareDifferenceVect_red);
                [valDiff_green idx_green_1] = min(SquareDifferenceVect_green);
                redSpot(frame, :) = [centroids_red(idx_red_1, 1), centroids_red(idx_red_1, 2)];
                greenSpot(frame, :) = [centroids_green(idx_green_1, 1), centroids_green(idx_green_1, 2)];
            end

        elseif (frame == first_frame) && ((no_red_led == 1) || (no_green_led ==1))
            first_frame = first_frame + 1;
        else
            
            if (~isempty(centroids_red))
                SquareDifferenceVect_red = [(centroids_red(:,1) - redSpot(frame-1, 1)).^2 + (centroids_red(:,2) - redSpot(frame-1, 2)).^2];
                [valDiff_red idx_red] = min(SquareDifferenceVect_red);
                if (valDiff_red < difference_threshold_red)
                    redSpot(frame, :) = [centroids_red(idx_red, 1), centroids_red(idx_red, 2)];
                    SquareDifference_red = [(redSpot(frame,1) - redSpot(frame-1, 1)).^2 + (redSpot(frame,2) - redSpot(frame-1, 2)).^2];
                else
                    redSpot(frame, :) = redSpot(frame-1, :);
                    SquareDifference_red = [(redSpot(frame,1) - redSpot(frame-1, 1)).^2 + (redSpot(frame,2) - redSpot(frame-1, 2)).^2];
                end
            else
                redSpot(frame, :) = redSpot(frame-1, :);
                SquareDifferenceVect_red = [(redSpot(frame,1) - redSpot(frame-1, 1)).^2 + (redSpot(frame,2) - redSpot(frame-1, 2)).^2];
            end
            
            if (~isempty(centroids_green))
                SquareDifferenceVect_green = [(centroids_green(:,1) - greenSpot(frame-1, 1)).^2 + (centroids_green(:,2) - greenSpot(frame-1, 2)).^2];
                [valDiff_green idx_green] = min(SquareDifferenceVect_green);
                greenSpot(frame, :) = [centroids_green(idx_green, 1), centroids_green(idx_green, 2)];
                if (valDiff_green < difference_threshold_green)
                    greenSpot(frame, :) = [centroids_green(idx_green, 1), centroids_green(idx_green, 2)];
                    SquareDifference_green = [(greenSpot(frame,1) - greenSpot(frame-1, 1)).^2 + (greenSpot(frame,2) - greenSpot(frame-1, 2)).^2];
                else
                    greenSpot(frame, :) = greenSpot(frame-1, :);
                    SquareDifference_green = [(greenSpot(frame,1) - greenSpot(frame-1, 1)).^2 + (greenSpot(frame,2) - greenSpot(frame-1, 2)).^2];
                end
            else
                greenSpot(frame, :) = greenSpot(frame-1, :);
                SquareDifference_green = [(greenSpot(frame,1) - greenSpot(frame-1, 1)).^2 + (greenSpot(frame,2) - greenSpot(frame-1, 2)).^2];
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
                [x_red, y_red] = getpts
                display('Choose green spot (use right-click)')
                [x_green, y_green] = getpts
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
        
        frame
            
    frame = frame + 1
        
    end
    
    HD_180_deg_total = [HD_180_deg_total HD_180_deg];

    
end

%         else
%             % Add case where either or both spots do not appear. Affect value 450 to frame. Save frame index here and affect value to HD after filtering
%             num_error_frames = num_error_frames + 1;
%             error_frame(num_error_frames) = frame;



% for i = 1:10
%     a(i) = load(horzcat('HD_180_deg_', num2str(i), '.mat'));
% end
% HD_180_deg = [];
% for i = 1:10
%     HD_180_deg = [HD_180_deg a(i).HD_180_deg];
% end





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