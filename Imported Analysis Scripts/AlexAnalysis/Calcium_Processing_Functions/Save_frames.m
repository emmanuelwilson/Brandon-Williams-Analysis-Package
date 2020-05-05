% clc
% close all
% clear all

% Open the VT1.avi demo movie that ships with MATLAB.
% First get the folder that it lives in.
folder = [];
movieFullFileName = [];
folder = fileparts(which(['behavCam1.avi'])); % Determine where video folder is.
movieFullFileName = fullfile(folder, ['behavCam1.avi']);
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



% Loop through the movie.
for frame = 1 : numberOfFrames
    % Extract the frame from the movie structure.
    thisFrame = read(videoObject, frame);
    
    %cla
    %imshow(thisFrame);
    imshow(thisFrame, 'Border', 'tight');
    
    filename = [sprintf('%05d',frame) '.tiff'];    %to save video
    saveas(gcf, filename);
    
    frame
    
end





%% Save video
% folder = fileparts(which(strcat('behavCam', num2str(11), '.avi')));  % msCam8   //// extraction videos 1,2 and 3 are msCam 3, 4 and 5 (3 dec)
% workingDir = folder;
%     imageNames = dir(fullfile(workingDir,['images_',num2str(1)],'*.tiff'));
% imageNames = {imageNames.name}';
%
% outputVideo = VideoWriter(fullfile(workingDir, strcat('extraction_test_', num2str(1), '.avi')));
% open(outputVideo)
%
% for ii = 1:length(imageNames)
%    img = imread(fullfile(workingDir,['images_',num2str(1)],imageNames{ii}));
%    writeVideo(outputVideo,img)
% end
%
% close(outputVideo)