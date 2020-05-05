%% Execute this before plot_head_direction_Zaki


% find dat files
datFiles = dir([pwd '\*.dat']);

%read timestamp information
for i=1:length(datFiles)
    if strcmp(datFiles(i).name,'timestamp.dat')
        fileID = fopen([pwd '\' datFiles(i).name],'r');
        dataArray = textscan(fileID, '%f%f%f%f%[^\n\r]', 'Delimiter', '\t', 'EmptyValue' ,NaN,'HeaderLines' ,1, 'ReturnOnError', false);
        camNum = dataArray{:, 1};
        frameNum = dataArray{:, 2};
        sysClock = dataArray{:, 3};
        buffer1 = dataArray{:, 4};
        clearvars dataArray;
        fclose(fileID);
        for j=max(camNum):-1:0   %% was for j=0:max(camNum) %%%%%%%%   Zaki   %%%%%%%%%%%%%% Change depending on the numbers of the behavioural and miniscope cameras
            %                 (frameNum(find(camNum==j,1,'last')) == ms.numFrames)
            %                 (sum(camNum==j) == ms.numFrames)
            if j == 1
                msSync.camNumber = j;
                msSync.time = sysClock(camNum == j);
                msSync.time(1) = 0;
                msSync.maxBufferUsed = max(buffer1(camNum==j));
                
            else
                behavSync.camNumber = j;
                behavSync.time = sysClock(camNum == j);
                behavSync.time(1) = 0;
                behavSync.maxBufferUsed = max(buffer1(camNum==j));
            end
        end
    end
end

frameMap = [];
halfWindow = 500;
if length(msSync.time) < length(behavSync.time)
    minRef = msSync.time;
    maxRef = behavSync.time;
else
    minRef = behavSync.time;
    maxRef = msSync.time;
end
for i = 1 : (length(minRef))
    timeDiff = [];
    for j = 1 : (halfWindow*2)+1
        if ((i + j - (halfWindow+1)) < 1)
            timeDiff(j) = NaN;
        elseif ((i + j - (halfWindow+1)) > length(maxRef))
            timeDiff(j) = NaN;
        else
            timeDiff(j) = abs(minRef(i) - maxRef(i + j - (halfWindow+1)));
   %         timeDiff(isnan(timeDiff(1,:)))=[];
            [minTimeDiff, indMinTimeDiff] = min(timeDiff);
        end
    end
        frameMap(i) = i + (indMinTimeDiff - (halfWindow+1));
end
frameMap = frameMap';
a = [frameMap msSync.time(1:length(frameMap)) behavSync.time(1:length(frameMap))];
