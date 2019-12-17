function ms = msGenerateVideoObj_TIFF(dirName, filePrefix)
%Organizes and spits out info and video files corresponding to "filePrefix" for further processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%   This function will read through all files inside the folder/directory
%   marked as "dirName". It will find all wanted .avi videofiles with a
%   specific name defined as "filePrefix", organize them via time stamp,
%   create a list of video objects of each .avi file and record frame
%   number, frame index and file # index. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Modified by Émmanuel Wilson

     v = VideoWriter('newfile.avi','Uncompressed AVI');
     open(v);
     for k=1:10      % assumes 10 images to write to file
         writeVideo(v,randi(255,100,200,'uint8'));
     end
     close(v);


    MAXFRAMESPERFILE = 1000; %This is set in the miniscope control software
    
    % find tif and dat files
    tifFiles = dir([dirName '\*.tif']);
    folder = dir(dirName);
    
    ms.numFiles = 0;        %Number of relevant .tif files in the folder
    ms.numFrames = 0;       %Number of frames within said videos
    ms.vidNum = [];         %Video index
    ms.frameNum = [];       %Frame number index
    ms.maxFramesPerFile = MAXFRAMESPERFILE; %finds the maximum number of frames contained in a single throughout all videos
    ms.dffframe = [];
    
    %find the total number of relevant video files
    for i=1:length(tifFiles)
        endIndex = strfind(tifFiles(i).name,'.tif');        %find the name of current .avi file
        if (~isempty(strfind(tifFiles(i).name,filePrefix)))
            ms.numFiles = max([ms.numFiles str2double(tifFiles(i).name((length(filePrefix)+1):endIndex))]);     % +1 count for relevant .avi files
        end
    end
    

    for i=1:ms.numFiles
        ms.vidObj{i} = VideoReader(folder(i).name);                     %Read .avi video file
        ms.vidNum = [ms.vidNum i*ones(1,ms.vidObj{i}.NumberOfFrames)];  %Store video index into ms for future use outside this fn
        ms.frameNum = [ms.frameNum 1:ms.vidObj{i}.NumberOfFrames];      %Current frame # in total
        ms.numFrames = ms.numFrames + ms.vidObj{i}.NumberOfFrames;      %Total number of frames
    end
    ms.height = ms.vidObj{1}.Height;        %video dimentions
    ms.width = ms.vidObj{1}.Width;
    
    camNum = [];
    frameNum = [];
    sysClock = [];
    buffer1 = [];
    frameTot1 = 0;
    frameTot0 = 0;

    idx = strfind(dirName, '_');
    idx2 = strfind(dirName,'\');
    if (length(idx) >= 4)
        ms.dateNum = datenum(str2double(dirName((idx(end-2)+1):(idx2(end)-1))), ... %year
            str2double(dirName((idx2(end-1)+1):(idx(end-3)-1))), ... %month
            str2double(dirName((idx(end-3)+1):(idx(end-2)-1))), ... %day
            str2double(dirName((idx2(end)+2):(idx(end-1)-1))), ...%hour
            str2double(dirName((idx(end-1)+2):(idx(end)-1))), ...%minute
            str2double(dirName((idx(end)+2):end)));%second
    end
end