fileID = fopen([pwd '\timestamp.dat'],'r');         %access dat files
dataArray = textscan(fileID, '%f%f%f%f%[^\n\r]', 'Delimiter', '\t', 'EmptyValue' ,NaN,'HeaderLines' ,1, 'ReturnOnError', false);    %read file and make sure it is not empty
camNum = dataArray{:, 1};       %camera number
frameNum = dataArray{:, 2};     %frame number
sysClock = dataArray{:, 3};     %system clock
buffer1 = dataArray{:, 4};      %buffer
clearvars dataArray;            %clear variables from dataArray
fclose(fileID);
j= 1;   %% was for j=0:max(camNum) %%%%%%%%   Zaki   %%%%%%%%%%%%%% Change depending on the numbers of the behavioural and miniscope cameras
%                 (frameNum(find(camNum==j,1,'last')) == ms.numFrames)
%                 (sum(camNum==j) == ms.numFrames)
if (sum(camNum==j)~=0)
camNumber = j;
time = sysClock(camNum == j);
time(1) = 0;
maxBufferUsed = max(buffer1(camNum==j));
msSync.camNumber = j;
msSync.time = sysClock(camNum == j);
msSync.time(1) = 0;
msSync.maxBufferUsed = max(buffer1(camNum==j));
end
j= 0;   %% was for j=0:max(camNum) %%%%%%%%   Zaki   %%%%%%%%%%%%%% Change depending on the numbers of the behavioural and miniscope cameras
%                 (frameNum(find(camNum==j,1,'last')) == ms.numFrames)
%                 (sum(camNum==j) == ms.numFrames)
if (sum(camNum==j)~=0)
camNumber = j;
time = sysClock(camNum == j);
time(1) = 0;
maxBufferUsed = max(buffer1(camNum==j));
msSync2.camNumber = j;
msSync2.time = sysClock(camNum == j);
msSync2.time(1) = 0;
msSync2.maxBufferUsed = max(buffer1(camNum==j));
end
timediff1 = diff(msSync.time);
timediff2 = diff(msSync2.time);
figure
plot(timediff2)
title('Cam0')
figure
plot(timediff1)
title('Cam1')