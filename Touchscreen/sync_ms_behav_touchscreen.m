%%Synchronize Miniscope with touchscreen data.

MAXFRAMESPERFILE = 1000; %This is set in the miniscope control software

%---------------------GRABBING DATA FROM CURRENT FOLDER--------------------
%---------------------(MINISCOPE AND TOUCH SCREEN FILES)-------------------
% find avi and dat files
aviFiles = dir([pwd '\*.avi']);
datFiles = dir([pwd '\*.dat']);
folder = dir(pwd);
filePrefix = 'msCam';
Prefix_TUNL = '20JAN_092804.dat'; 

numFiles = 0;        %Number of relevant .avi files in the folder
numFrames = 0;       %Number of frames within said videos
vidNum = [];         %Video index
frameNum = [];       %Frame number index
maxFramesPerFile = MAXFRAMESPERFILE; %finds the maximum number of frames contained in a single throughout all videos
dffframe = [];

%find the total number of relevant video files
for i=1:length(aviFiles)
    endIndex = strfind(aviFiles(i).name,'.avi');        %find the name of current .avi file
    if (~isempty(strfind(aviFiles(i).name,filePrefix)))
        numFiles = max([numFiles str2double(aviFiles(i).name((length(filePrefix)+1):endIndex))]);     % +1 count for relevant .avi files
        vidObj{i} = VideoReader(aviFiles(i).name);                     %Read .avi video file
        vidNum = [vidNum i*ones(1,vidObj{i}.NumberOfFrames)];  %Store video index into ms for future use outside this fn
        frameNum = [frameNum 1:vidObj{i}.NumberOfFrames];      %Current frame # in total
        numFrames = numFrames + vidObj{i}.NumberOfFrames;      %Total number of frames
    end
end

% height = vidObj{1}.Height;        %video dimentions
% width = vidObj{1}.Width;

%read timestamp information
for i=1:length(datFiles)
    if strcmp(datFiles(i).name,'timestamp.dat')
        fileID = fopen([pwd '\' datFiles(i).name],'r');         %access dat files
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
                if ((frameNum(find(camNum==j,1,'last')) == numFrames) && (sum(camNum==j) == numFrames))
                    camNumber = j;
                    time = sysClock(camNum == j);
                    time(1) = 0;
                    maxBufferUsed = max(buffer1(camNum==j));
                    msSync.camNumber = j;
                    msSync.time = sysClock(camNum == j);
                    msSync.time(1) = 0;
                    msSync.maxBufferUsed = max(buffer1(camNum==j));
                else
                    display(['Problem matching up timestamps for ' pwd]);
                    behavSync.camNumber = j;
                    behavSync.time = sysClock(camNum == j);
                    behavSync.time(1) = 0;
                    behavSync.maxBufferUsed = max(buffer1(camNum==j));
                end
            end
        end
    elseif strcmp(datFiles(i).name,Prefix_TUNL) %currently you need to put the name of the dat file here
        fileID = fopen([pwd '\' datFiles(i).name],'r');         %access dat files
        datArray = textscan(fileID,  '%q%q%q%q%q%q%q%[^\n\r]', 'Delimiter', '\t', 'EmptyValue' ,NaN,'HeaderLines' ,1, 'ReturnOnError', false);    %read file and make sure it is not empty
        Evnt_Time = datArray{:, 1}; %time is in: ___
        Evnt_Name = datArray{:,3};
        Item_Name = datArray{:, 4};
        Alias_Name = datArray{:,5}; %What is alias name?
        Group_ID = datArray{:, 6};
        Arg1 = datArray{:,8};        
        clearvars datArray;            %clear variables from dataArray
        fclose(fileID);
    end
end

%-----------------------CONVERTING TOUCH SCREEN VARIABLES------------------
eventTime = zeros(size(Evnt_Time));   %Touchscreen event timestamps
eventname = strings(size(eventTime));
item_name = strings(size(eventTime));  %Touchscreen event name
aliasname = eventname;
groupID = eventTime;  %Touchscreen event group ID
arg = eventname;

%convert cell arrays into arrays
for t = 1 : length(eventTime)                                               
    eventTime(t) = str2double(Evnt_Time{t});
    eventname(t) = Evnt_Name{t};
    item_name(t) = Item_Name{t};
    aliasname(t) = Alias_Name{t};
    groupID(t) = str2double(Group_ID{t});
    arg(t) = Arg1{t};    
end
%--------------------------------------------------------------------------

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
        elseif ~isempty(find(contains(eventname(group),'Touch Down Event'),1)) && ~isempty(find(contains(item_name(group),'Incorrect'),1)) && ~isempty(find(contains(item_name(group),'Start Delay'),1))
            events(eventInd(i),2) = timeMap(eventInd(i));
            events(eventInd(i),3) = 'Nose-Poke Incorrect';
            [a,b]=strtok(arg(group(find(contains(eventname(group),'Touch Down Event')+i-1,1))),['Image' 'Position']);            
            events(eventInd(i),4) = strtok(a);
            events(eventInd(i),6) = 1;
        elseif ~isempty(find(contains(eventname(group),'Touch Down Event'),1)) && ~isempty(find(contains(item_name(group),'Correct'),1)) && ~isempty(find(contains(item_name(group),'Start Delay'),1))
            events(eventInd(i),2) = timeMap(eventInd(i));
            events(eventInd(i),3) = 'Nose-Poke Correct';
            [a,b]=strtok(arg(group(find(contains(eventname(group),'Touch Down Event')+i-1,1))),['Image' 'Position']);            
            events(eventInd(i),4) = strtok(a);
            events(eventInd(i),6) = 1;
        elseif ~isempty(find(contains(eventname(group),'Touch Down Event'),1)) && ~isempty(find(contains(item_name(group),'Incorrect'),1))
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
            group2 = find(eventInd == eventInd(max(group)+1));
            if ~isempty(find(contains(item_name(group2),'Incorrect'),1)) && ~isempty(find(contains(item_name(group),'Start Delay'),1))
                events(eventInd(i),2) = timeMap(eventInd(i));
                events(eventInd(i),3) = 'Nose-Poke Incorrect';
                [a,b]=strtok(arg(group(find(contains(eventname(group),'Touch Down Event')+i-1,1))),['Image' 'Position']);            
                events(eventInd(i),4) = strtok(a);
                events(eventInd(i),6) = 1;
            elseif ~isempty(find(contains(item_name(group2),'Correct'),1))&& ~isempty(find(contains(item_name(group),'Start Delay'),1))
                events(eventInd(i),2) = timeMap(eventInd(i));
                events(eventInd(i),3) = 'Nose-Poke Correct';
                [a,b]=strtok(arg(group(find(contains(eventname(group),'Touch Down Event')+i-1,1))),['Image' 'Position']);            
                events(eventInd(i),4) = strtok(a);
                events(eventInd(i),6) = 1;
            elseif ~isempty(find(contains(item_name(group2),'Incorrect'),1)) && isempty(find(contains(item_name(group2),'Next trial'),1)) 
                events(eventInd(i),2) = timeMap(eventInd(i));
                events(eventInd(i),3) = 'Nose-Poke Incorrect';
                [a,b]=strtok(arg(group(find(contains(eventname(group),'Touch Down Event')+i-1,1))),['Image' 'Position']);            
                events(eventInd(i),4) = strtok(a);
            elseif ~isempty(find(contains(item_name(group2),'Correct'),1)) && isempty(find(contains(item_name(group2),'Next trial'),1)) 
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

save('msTouchSync.mat', 'events', 'frameMap','eventTime','eventInd', 'item_name','timeMap','groupID','arg');