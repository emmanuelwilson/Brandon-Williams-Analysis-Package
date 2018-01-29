function ms = msGenerateVideoObj_Emm(dirName, filePrefix)
%Organizes and spits out info and video files corresponding to "filePrefix" for further processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%   This function will read through all files inside the folder/directory
%   marked as "dirName". It will find all wanted .avi videofiles with a
%   specific name defined as "filePrefix", organize them via time stamp,
%   create a list of video objects of each .avi file and record frame
%   number, frame index and file # index. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Modified by Émmanuel Wilson

    MAXFRAMESPERFILE = 1000; %This is set in the miniscope control software
    
    % find avi and dat files
    aviFiles = dir([dirName '\*.avi']);
    datFiles = dir([dirName '\*.dat']);
    folder = dir(dirName);
    
    ms.numFiles = 0;        %Number of relevant .avi files in the folder
    ms.numFrames = 0;       %Number of frames within said videos
    ms.vidNum = [];         %Video index
    ms.frameNum = [];       %Frame number index
    ms.maxFramesPerFile = MAXFRAMESPERFILE; %finds the maximum number of frames contained in a single throughout all videos
    ms.dffframe = [];
    
    %find the total number of relevant video files
    for i=1:length(aviFiles)
        endIndex = strfind(aviFiles(i).name,'.avi');        %find the name of current .avi file
        if (~isempty(strfind(aviFiles(i).name,filePrefix)))
            ms.numFiles = max([ms.numFiles str2double(aviFiles(i).name((length(filePrefix)+1):endIndex))]);     % +1 count for relevant .avi files
        end
    end
    
    Timestampmin = 1000000;     %Arbitrary VERY large number for min value storage, must be larger than datenum stamp in folder
    Timestamp = 0;              %Total number of frames observed (including previous video itterations)
    o = NaN(1,ms.numFiles);     %index movie file order
    anomilynames = NaN(1,4);
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
    if folder(o(1,length(o(1,:)))).bytes > folder(o(1,length(o(1,:))-1)).bytes
        temp = o(1,length(o(1,:)));
        o(1,length(o(1,:)))= o(1,length(o(1,:))-1);
        o(1,length(o(1,:))-1)= temp;
    end
    %generate a vidObj for each video file. Also calculate total frames
    for i=1:ms.numFiles
        j = o(1,i);                                                     %call on .avi files chronologically   
        ms.vidObj{i} = VideoReader(folder(j).name);                     %Read .avi video file
        ms.vidNum = [ms.vidNum i*ones(1,ms.vidObj{i}.NumberOfFrames)];  %Store video index into ms for future use outside this fn
        ms.frameNum = [ms.frameNum 1:ms.vidObj{i}.NumberOfFrames];      %Current frame # in total
        ms.numFrames = ms.numFrames + ms.vidObj{i}.NumberOfFrames;      %Total number of frames
    end
    ms.height = ms.vidObj{1}.Height;        %video dimentions
    ms.width = ms.vidObj{1}.Width;
    
    %read timestamp information
    for i=1:length(datFiles)
        if strcmp(datFiles(i).name,'timestamp.dat')
            fileID = fopen([dirName '\' datFiles(i).name],'r');         %access dat files
            dataArray = textscan(fileID, '%f%f%f%f%[^\n\r]', 'Delimiter', '\t', 'EmptyValue' ,NaN,'HeaderLines' ,1, 'ReturnOnError', false);    %read file and make sure it is not empty
            camNum = dataArray{:, 1};       %camera number
            frameNum = dataArray{:, 2};     %frame number 
            sysClock = dataArray{:, 3};     %system clock
            buffer1 = dataArray{:, 4};      %buffer
            clearvars dataArray;            %clear variables from dataArray
            fclose(fileID);
            for j=max(camNum):-1:0   %% was for j=0:max(camNum) %%%%%%%%   Zaki   %%%%%%%%%%%%%% Change depending on the numbers of the behavioural and miniscope cameras 
%                 (frameNum(find(camNum==j,1,'last')) == ms.numFrames)
%                 (sum(camNum==j) == ms.numFrames)
                if (sum(camNum==j)~=0)
                if ((frameNum(find(camNum==j,1,'last')) == ms.numFrames) && (sum(camNum==j) == ms.numFrames))
                    ms.camNumber = j;                   
                    ms.time = sysClock(camNum == j);
                    ms.time(1) = 0;
                    ms.maxBufferUsed = max(buffer1(camNum==j));
                
                else
                    display(['Problem matching up timestamps for ' dirName]);        
                end
                end
            end
        end
%         
%         %read in and store animal name
%         if strcmp(datFiles(i).name, 'settings_and_notes.dat')
%             fileID = fopen([dirName '\' datFiles(i).name],'r');
%             textscan(fileID, '%[^\n\r]', 1, 'ReturnOnError', false);
%             dataArray = textscan(fileID, '%s%s%s%s%[^\n\r]', 1, 'Delimiter', '\t', 'ReturnOnError', false);
%             ms.animal = dataArray(:,1);
%             ms.animal = ms.animal{1};
%         end
    end
    
%     
%     %figure out date and time of recording if that information if available
%     %in folder path
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