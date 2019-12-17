
dirName = pwd;
    datFiles = dir([dirName '\*.dat']);
    for i=1:length(datFiles)
        if strfind(datFiles(i).name,'timestamp')
            fileID = fopen([dirName '\' datFiles(i).name],'r');         %access dat files
            dataArray = textscan(fileID, '%f%f%f%f%[^\n\r]', 'Delimiter', '\t', 'EmptyValue' ,NaN,'HeaderLines' ,1, 'ReturnOnError', false);    %read file and make sure it is not empty
            camNumtemp = dataArray{:, 1};       %camera number
            frameNumtemp = dataArray{:, 2};     %frame number 
            sysClocktemp = dataArray{:, 3};     %system clock
            buffer1temp = dataArray{:, 4};      %buffer
            clearvars dataArray;            %clear variables from dataArray
            fclose(fileID);
            if i == 1
                camNum = camNumtemp ;
                frameNum =  frameNumtemp;
                sysClock = sysClocktemp;
                buffer1 = buffer1temp;
            else
                camNum = [camNum; camNumtemp] ;
                frameNum = [frameNum; frameNumtemp];
                sysClock = [sysClock; sysClock(length(sysClock))+sysClocktemp];
                buffer1 = [buffer1; buffer1temp];
            end
            %shouldn't we take max and not length because of repeating
            %frame numbers ? 
            if mod(i,2)
                frameTot1 = frameTot1 + frameNum(max(frameNum));
                sysClockTot = [sysClockTot; max(sysClockTot(2:end)) + sysClock];
                frameTot0 = frameTot0 + frameNum(max(frameNum)-1);
            else
                frameTot0 = frameTot0 + frameNum(max(frameNum));
                sysClockTot = [sysClockTot; max(sysClockTot(3:end)) + sysClock];
                frameTot1 = frameTot1 + frameNum(max(frameNum)-1);
            end
        end

    end
    for j=max(camNum):-1:0   %% was for j=0:max(camNum) %%%%%%%%   Zaki   %%%%%%%%%%%%%% Change depending on the numbers of the behavioural and miniscope cameras
        %                 (frameNum(find(camNum==j,1,'last')) == ms.numFrames)
        %                 (sum(camNum==j) == ms.numFrames)
        if (sum(camNum==j)~=0)
            if (frameTot1 == ms.numFrames) && (sum(camNum==j) == ms.numFrames)
                ms.camNumber = j;
                ms.time = sysClock(camNum == j);
                ms.time(1) = 0;
                ms.maxBufferUsed = max(buffer1(camNum==j));
            elseif (frameTot0 == ms.numFrames) && (sum(camNum==j) == ms.numFrames)
                ms.camNumber = j;
                ms.time = sysClock(camNum == j);
                ms.time(1) = 0;
                ms.maxBufferUsed = max(buffer1(camNum==j));
            else
                display(['Problem matching up timestamps for ' dirName]);
            end
        end
    end