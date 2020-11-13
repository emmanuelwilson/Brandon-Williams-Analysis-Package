%%Synchronize Miniscope with touchscreen data.
function msTouchSync = sync_ms_behav_touchscreen_V5(folderpath)

MAXFRAMESPERFILE = 1000; %This is set in the miniscope control software

%---------------------GRABBING DATA FROM CURRENT FOLDER--------------------
%---------------------(MINISCOPE AND TOUCH SCREEN FILES)-------------------
% find avi and dat files
aviFiles = dir([folderpath '\*.avi']);
csvFiles = dir([folderpath '\*.csv']);

aviFilesbehav = dir([folderpath '\..\BehavCam*\*.avi']);
csvFilesbehav = dir([folderpath '\..\BehavCam*\*.csv']);

filePrefix = '';
%Prefix_TUNL = 'CI28759-3_TUNL_J20_CA1_TUNL Mouse Exp 1 Stage 1 S3TTL_D4sec_232.csv';

numFrames = 0;       %Number of frames within said videos
numFramesbehav = 0;       %Number of frames within said videos

%find the total number of relevant video files
for i=1:length(aviFiles)    
    try
        vidObj{i} = VideoReader([aviFiles(i).folder '\' aviFiles(i).name]);                     %Read .avi video file
        numFrames = numFrames + vidObj{i}.NumberOfFrames;      %Total number of frames
    catch
        fprintf(['Calcium Video ' aviFiles(i).name ' is Invalid'])
        stop
    end
end
vidObj = [];

%find the total number of relevant Behav video files
for i=1:length(aviFilesbehav)    
    try
    vidObj{i} = VideoReader([aviFilesbehav(i).folder '\' aviFilesbehav(i).name]);                     %Read .avi video file
    numFramesbehav = numFramesbehav + vidObj{i}.NumberOfFrames;      %Total number of frames
    catch
        fprintf(['Behaviour Video ' aviFilesbehav(i).name ' is Invalid'])
        stop
    end
end

%read timestamp information
for i=1:length(csvFiles)
    if strcmp(csvFiles(i).name,'timeStamps.csv')
        fileID = fopen([csvFiles(i).folder '\' csvFiles(i).name],'r');         %access dat files
        dataArray = textscan(fileID, '%f%f%f%f%[^\n\r]', 'Delimiter', ',', 'EmptyValue' ,NaN,'HeaderLines' ,1, 'ReturnOnError', false);    %read file and make sure it is not empty        
        frameNum = dataArray{:, 1};     %frame number
        sysClock = dataArray{:, 2};     %system clock
        buffer1 = dataArray{:, 3};      %buffer
        clearvars dataArray;            %clear variables from dataArray
        fclose(fileID);                           
        if ((frameNum(end))+1 == numFrames) && (length(frameNum) == numFrames)            
            time = sysClock;
            time(1) = 0;                       
            msSync.time = sysClock;
            msSync.time(1) = 0;
            msSync.maxBufferUsed = max(buffer1);
        else
            fprintf('Calcium video and Timestamp frame number missmatch');
        end
    end
end
for i=1:length(csvFilesbehav)
    if strcmp(csvFilesbehav(i).name,'timeStamps.csv')
        fileID = fopen([csvFilesbehav(i).folder '\' csvFilesbehav(i).name],'r');         %access dat files
        dataArray = textscan(fileID, '%f%f%f%f%[^\n\r]', 'Delimiter', ',', 'EmptyValue' ,NaN,'HeaderLines' ,1, 'ReturnOnError', false);    %read file and make sure it is not empty        
        frameNum = dataArray{:, 1};     %frame number
        sysClock = dataArray{:, 2};     %system clock
        buffer1 = dataArray{:, 3};      %buffer
        clearvars dataArray;            %clear variables from dataArray
        fclose(fileID);                           
        if ((frameNum(end))+1 == numFramesbehav) && (length(frameNum) == numFramesbehav)            
            time = sysClock;
            time(1) = 0;
            behavSync.time = sysClock;
            behavSync.time(1) = 0;
            behavSync.maxBufferUsed = max(buffer1);
        else
            fprintf('Behaviour video and Timestamp frame number missmatch');
        end
    end
end

i = 1;
Tschedule = readtable([folderpath '\' csvFiles(i).name]);
s = size(Tschedule);
while s(2) ~= 17    
    if i > length(csvFiles) || s(2) == 17
        break;
    else
        i = i +1;
        Tschedule = readcell([folderpath '\' csvFiles(i).name]);
        s = size(Tschedule);
    end        
end

% Tschedule(1:16,:) = [];
Evnt_Time = Tschedule(:, 1); %time is in: ___
Evnt_Name = Tschedule(:, 3);
Item_Name = Tschedule(:, 4);
Alias_Name = Tschedule(:, 5); %What is alias name?
Group_ID = Tschedule(:, 6);
Arg1 = Tschedule(:,9);

%-----------------------CONVERTING TOUCH SCREEN VARIABLES------------------
eventTime = zeros(size(Evnt_Time));   %Touchscreen event timestamps
eventname = strings(size(eventTime));
item_name = strings(size(eventTime));  %Touchscreen event name
aliasname = eventname;
groupID = eventTime;  %Touchscreen event group ID
arg = eventname;


%convert cell arrays into arrays
for t = 1 : length(eventTime)                                               
    eventTime(t) = Evnt_Time{t,1};
    eventname(t) = Evnt_Name{t,1};
    item_name(t) = Item_Name{t,1};
    aliasname(t) = Alias_Name{t,1};
    groupID(t) = Group_ID{t,1};
    arg(t) = Arg1{t,1};    
end

frameMap = [];
timeMap = []; 
halfWindow = 500; %what's half window?? half 100 frames from miniscope videos?

%select shortest video to be baseline frame number index
%if behav is longer than miniscope, use miniscope time as minimum reference
%and behav as maximum reference?
if length(msSync.time) < length(behavSync.time)
    minRef = msSync.time; %what is this for? 
    maxRef = behavSync.time;    
    time = behavSync.time/1000; %why are you using behav for the time?
else %if time is longer, use miniscope time as reference? 
    minRef = behavSync.time;
    maxRef = msSync.time;   
    time = time/1000;
end

%filling the frame map
for i = 1 :(length(minRef))
    timeDiff = [];
    for j = 1 : (halfWindow*2)+1
        if ((i + j - (halfWindow+1)) < 1)
            timeDiff(j) = NaN;
        elseif ((i + j - (halfWindow+1)) > length(maxRef))
            timeDiff(j) = NaN;
        else
            timeDiff(j) = abs(minRef(i) - maxRef(i + j - (halfWindow+1)));
            [minTimeDiff, indMinTimeDiff] = min(timeDiff);
        end
    end
    frameMap(i) = i + (indMinTimeDiff - (halfWindow+1));
    timeMap(i) = time(i + (indMinTimeDiff - (halfWindow+1)));
end

frameMap = frameMap';
timeMap = timeMap';
a = [frameMap msSync.time(1:length(frameMap)) behavSync.time(1:length(frameMap))];

loc = strfind(item_name,'TTL');                                             %find the TTL initialization
startTime = NaN;        
startInd = NaN;
endTime = NaN;
endInd = NaN;
firstTTL = 0;
%find miniscope start point in behaviour
for t = 1 : length(eventTime)
    if ~isempty(loc{t}) && firstTTL ==0
        startTime = eventTime(t);
        startInd = t;
        firstTTL = 1;
    elseif ~isempty(loc{t}) && firstTTL >0
        endTime = eventTime(t);
        endInd = t;        
    end
    
end
eventInd = [];%nan(length(eventTime)-startInd,1);
timeMap(:,1)  = timeMap(:,1) + startTime;                                   %Add the miniscope start time to miniscope time
%behaviour event time to miniscope time 
for t = 1 : length(eventTime)
    if eventTime(t) <= max(timeMap)
        [minDiff, ind ] = min(abs(eventTime(t) - timeMap));
        eventInd(t) = ind;
    end
end
eventInd = eventInd';

events = strings(length(frameMap),6);
events(:,1) = string(frameMap);
events(eventInd(startInd),2) = timeMap(eventInd(startInd));
events(eventInd(startInd),3) = item_name(startInd);
ind = 1;

% Looking through the strings to find and save the wanted events
%events is split as follows: 
%events(:,1) = frameMap
%events(:,2) = timestamp of event
%events(:,3) = event name
%events(:,4) = nose-poke position
%events(:,5) = Trial Start/stop, indicated by a 1
%events(:,6) = Delay period start/stop (1 = start, 2 = stop)

for i = startInd : length(eventInd)
    if abs(ind - eventInd(i))>0
        group = find(eventInd == eventInd(i));
        [c,d] = strtok(item_name(group(find(contains(item_name(group),'Correction_Trial'),1))),'Correction_Trials');
        [a1,a2]=strtok(arg(group(find(contains(item_name(group),'Correction_Trial'),1))),'Value');
        if ~isempty(find(contains(item_name(group),'BIRBeam'),1)) && ~isempty(find(contains(item_name(group),'Correction_Trial'),1))&& c=='' && d=='' && contains(a1,'1') && ~isempty(find(contains(eventname(group),'Display Image'),1)) && ~isempty(find(contains(aliasname(group),'Training'),1))
            events(eventInd(i),2) = timeMap(eventInd(i));
            events(eventInd(i),3) = 'Display Image & BIRBeam#1 & correction Trial Start';
            [a,b]=strtok(arg(group(find(contains(eventname(group),'Display Image'),1))),['Image' 'Position']);
            events(eventInd(i),4) = strtok(a);
        elseif ~isempty(find(contains(eventname(group),'Display Image'),1)) && ~isempty(find(contains(aliasname(group),'Training'),1)) && ~isempty(find(contains(item_name(group),'BIRBeam'),1))
            events(eventInd(i),2) = timeMap(eventInd(i));
            events(eventInd(i),3) = 'BIRBeam & Display Image';
            [a,b]=strtok(arg(group(find(contains(eventname(group),'Display Image'),1))),['Image' 'Position']);
            events(eventInd(i),4) = strtok(a);
        elseif ~isempty(find(contains(eventname(group),'Display Image'),1)) && ~isempty(find(contains(aliasname(group),'Training'),1))
            events(eventInd(i),2) = timeMap(eventInd(i));
            events(eventInd(i),3) = 'Display Image';
            [a,b]=strtok(arg(group(find(contains(eventname(group),'Display Image'),1))),['Image' 'Position']);
            events(eventInd(i),4) = strtok(a);
        elseif ~isempty(find(contains(item_name(group),'BIRBeam'),1)) && ~isempty(find(contains(item_name(group),'Correction_Trial'),1)) && contains(a1,'1') && c=='' && d==''
            events(eventInd(i),2) = timeMap(eventInd(i));
            events(eventInd(i),3) = 'BIRBeam#1 & correction Trial Start';
        elseif ~isempty(find(contains(item_name(group),'BIRBeam'),1))
            events(eventInd(i),2) = timeMap(eventInd(i));
            events(eventInd(i),3) = 'BIRBeam#1';                
        elseif ~isempty(find(contains(eventname(group),'Touch Down Event'),1)) && ~isempty(find(strcmp(item_name(group),'Incorrect'),1)) && ~isempty(find(contains(item_name(group),'Start Delay'),1))
            events(eventInd(i),2) = timeMap(eventInd(i));
            events(eventInd(i),3) = 'Nose-Poke Incorrect';
            [a,b]=strtok(arg(group(find(contains(eventname(group),'Touch Down Event')+i-1,1))),['Image' 'Position']);            
            events(eventInd(i),4) = strtok(a);
            events(eventInd(i),6) = 1;
        elseif ~isempty(find(contains(eventname(group),'Touch Down Event'),1)) && ~isempty(find(strcmp(item_name(group),'Correct'),1)) && ~isempty(find(contains(item_name(group),'Start Delay'),1))
            events(eventInd(i),2) = timeMap(eventInd(i));
            events(eventInd(i),3) = 'Nose-Poke Correct';
            [a,b]=strtok(arg(group(find(contains(eventname(group),'Touch Down Event')+i-1,1))),['Image' 'Position']);            
            events(eventInd(i),4) = strtok(a);
            events(eventInd(i),6) = 1;
        elseif ~isempty(find(contains(eventname(group),'Touch Down Event'),1)) && ~isempty(find(strcmp(item_name(group),'Incorrect'),1)) || (i < length(eventInd) && ~isempty(find(contains(eventname(group),'Touch Down Event'),1)) && ~isempty(find(contains(item_name(find(eventInd == eventInd(i+1))),'Incorrect'),1)) && isempty(find(contains(eventname(find(eventInd == eventInd(i+1))),'Touch Down Event'),1)) )|| (i < length(eventInd)-1 && ~isempty(find(contains(eventname(group),'Touch Down Event'),1)) && isempty(find(contains(eventname(find(eventInd == eventInd(i+1))),'Touch Down Event'),1)) && ~isempty(find(contains(item_name(find(eventInd == eventInd(i+2))),'Incorrect'),1)) && isempty(find(contains(eventname(find(eventInd == eventInd(i+2))),'Touch Down Event'),1)))
            events(eventInd(i),2) = timeMap(eventInd(i));
            events(eventInd(i),3) = 'Nose-Poke Incorrect';
            [a,b]=strtok(arg(group(find(contains(eventname(group),'Touch Down Event')+i-1,1))),['Image' 'Position']);            
            events(eventInd(i),4) = strtok(a);
        elseif ~isempty(find(contains(eventname(group),'Touch Down Event'),1)) && ~isempty(find(contains(item_name(group),'Correct'),1))
            events(eventInd(i),2) = timeMap(eventInd(i));
            events(eventInd(i),3) = 'Nose-Poke Correct';
            [a,b]=strtok(arg(group(find(contains(eventname(group),'Touch Down Event')+i-1,1))),['Image' 'Position']);            
            events(eventInd(i),4) = strtok(a);
        elseif ~isempty(find(contains(eventname(group),'Touch Down Event'),1))
            if max(group) < length(eventInd)
                group2 = find(eventInd == eventInd(max(group)+1));
            else
                group2 = [];
            end
            if ~isempty(group2) && ~isempty(find(strcmp(item_name(group2),'Incorrect'),1)) && ~isempty(find(contains(item_name(group),'Start Delay'),1))
                events(eventInd(i),2) = timeMap(eventInd(i));
                events(eventInd(i),3) = 'Nose-Poke Incorrect';
                [a,b]=strtok(arg(group(find(contains(eventname(group),'Touch Down Event')+i-1,1))),['Image' 'Position']);            
                events(eventInd(i),4) = strtok(a);
                events(eventInd(i),6) = 1;
            elseif ~isempty(group2) && ~isempty(find(strcmp(item_name(group2),'Correct'),1))&& ~isempty(find(contains(item_name(group),'Start Delay'),1))
                events(eventInd(i),2) = timeMap(eventInd(i));
                events(eventInd(i),3) = 'Nose-Poke Correct';
                [a,b]=strtok(arg(group(find(contains(eventname(group),'Touch Down Event')+i-1,1))),['Image' 'Position']);            
                events(eventInd(i),4) = strtok(a);
                events(eventInd(i),6) = 1;
            elseif ~isempty(group2) && ~isempty(find(strcmp(item_name(group2),'Incorrect'),1)) && isempty(find(contains(item_name(group2),'Next trial'),1)) 
                events(eventInd(i),2) = timeMap(eventInd(i));
                events(eventInd(i),3) = 'Nose-Poke Incorrect';
                [a,b]=strtok(arg(group(find(contains(eventname(group),'Touch Down Event')+i-1,1))),['Image' 'Position']);            
                events(eventInd(i),4) = strtok(a);
            elseif ~isempty(group2) && ~isempty(find(strcmp(item_name(group2),'Correct'),1)) && isempty(find(contains(item_name(group2),'Next trial'),1)) 
                events(eventInd(i),2) = timeMap(eventInd(i));
                events(eventInd(i),3) = 'Nose-Poke Correct';
                [a,b]=strtok(arg(group(find(contains(eventname(group),'Touch Down Event')+i-1,1))),['Image' 'Position']);            
                events(eventInd(i),4) = strtok(a);
            elseif ~isempty(find(contains(item_name(group),'Start Delay'),1))           
                events(eventInd(i),2) = timeMap(eventInd(i));
                events(eventInd(i),3) = 'Nose-Poke';
                [a,b]=strtok(arg(group(find(contains(eventname(group),'Touch Down Event')+i-1,1))),['Image' 'Position']);
                events(eventInd(i),4) = strtok(a);
                events(eventInd(i),6) = 1;
            else
                if timeMap(eventInd(i)) > 150
                end
                events(eventInd(i),2) = timeMap(eventInd(i));
                events(eventInd(i),3) = 'Nose-Poke';
                [a,b]=strtok(arg(group(find(contains(eventname(group),'Touch Down Event')+i-1,1))),['Image' 'Position']);
                events(eventInd(i),4) = strtok(a);
            end        
        elseif ~isempty(find(contains(eventname(group),'Touch Up'),1))
            events(eventInd(i),2) = timeMap(eventInd(i));
            events(eventInd(i),3) = eventname(i);
            [a,b]=strtok(arg(group(find(contains(eventname(group),'Touch Up Event')+i-1,1))),['Image' 'Position']);            
            events(eventInd(i),4) = strtok(a);        
        else
            events(eventInd(i),2) = timeMap(eventInd(i));
            events(eventInd(i),3) = item_name(i);
        end
        if ~isempty(find(contains(item_name(group),'Correction_Trial'),1))&& contains(a1,'1') && c=='' && d==''
            events(eventInd(i),2) = timeMap(eventInd(i));
            events(eventInd(i),3) = 'correction Trial';
        elseif ~isempty(find(contains(item_name(group),'Start Delay'),1))
            events(eventInd(i),2) = timeMap(eventInd(i));
            events(eventInd(i),6) = 1;
        elseif ~isempty(find(contains(item_name(group),'Delay End'),1)) 
            events(eventInd(i),2) = timeMap(eventInd(i));
            events(eventInd(i),6) = 2;
        end
        if ~isempty(find(contains(item_name(group),'Next trial'),1)) &&   strcmp(item_name(group(find(contains(item_name(group),'Next trial'),1))), 'Next trial')          
            events(eventInd(i),5) = '1';
        end
    end
end

msTouchSync = struct;
msTouchSync.events = events;
msTouchSync.frameMap = frameMap;
msTouchSync.eventTime = eventTime;
msTouchSync.eventInd = eventInd;
msTouchSync.item_name = item_name;
msTouchSync.timeMap = timeMap;
msTouchSync.groupID = groupID;
msTouchSync.arg = arg;

msTouchSync.folder = folderpath;

save('msTouchSync.mat', 'msTouchSync');
end
