function ms = msSelectROIs(ms)
%This funtion delivers the user with a frame in a video segment, prompts
%the user to drag the mouse and select a rectangular region. After which
%asks if an additional region is to be selected. The selected region is
%saved into rect and later into AlignmentROI to be used in functions
%msAlignmentFFt etc. 

    numROIs=0;                                                      %Here refFrameNumber takes in the number
    userInput = 'Y';                                                %of frames in the ms structure and finds a middle frame number. 
    refFrameNumber = ceil(ms.numFrames/2);                          %example, 600 frames and we use the 300th frame. We then find 
    refFrame = msReadFrame(ms,refFrameNumber,false,false,false);    %that numbered frame in a video found in msReadFrame and set it to refFrame. 
                                                                    %i.e. we take the 300th from a selected video and assign it to frame in refFrame.     
    
    % starts by checking if alignmentROIs exists in ms, if true statment proceeds 
    % and displays the previously found frame in                                                                
    % specified fluroescent colours in ms.
    % in the for loop, as long as ROINum is 1:the size of the second number
    % in alignmentROI the loop will make a rectangle with set settings.  
    % Once the loop has finshed userInput asks for the user of the program
    % to reset the session or not 
    if (isfield(ms,'alignmentROI')) 
        imshow(uint8(refFrame), [min(ms.minFluorescence) max(ms.maxFluorescence)])
        hold on
        for ROINum = 1:size(ms.alignmentROI,2)
            rectangle('Position', ms.alignmentROI(:,ROINum),'LineWidth',2,'LineStyle','--');
        end
        userInput = upper(input('Session already has alignment ROIs. Reset ROIs? (Y/N)','s'));
    end

    %if the user of the program enters Y for yes then the session continues
    %through the statment. Else if the user had entered N, session would
    %have been terminated. 
    
    %once temp is defined, idx identifies if it is within ms, if yes idx =1
    %and the next line, ms removes temp from within its structure. 
    if strcmp(userInput,'Y')
        ms.alignmentROI = [];
        temp = {'hShift','wShift','alignedWidth','alignedHeight'};
        idx = isfield(ms,temp);
        ms = rmfield(ms,temp(idx));

        imshow(uint8(refFrame))     %converts refFrame to unit8 and displays the frame 
        hold on                     %this will display with the previous frame
        
        %While the user continues to input a Y (yes) into the program, the
        %loop will continue. 
        while (strcmp(userInput,'Y'))
            numROIs = numROIs+1; %produces a counter for number of ROIs 
            display(['Select ROI #' num2str(numROIs)])  
            rect = getrect();%allows the user to select a rectangle in the current display 
                             %and sets rect to that selected space 
            rect(3) = rect(3) - mod(rect(3),2); 
            rect(4) = rect(4) - mod(rect(4),2);

            %now we assign the values placed in rect to alignmentROI for
            %use in msAlignmentFTT 
            ms.alignmentROI(:,numROIs) = rect; %uint16([rect(1) rect(1)+rect(3) rect(2) rect(2)+rect(4)]);
            rectangle('Position',rect,'LineWidth',2); %changing rectangle to accomodate alignment to rect but adding rect in its place 
            userInput = upper(input('Select another ROI? (Y/N)','s')); %promptes the user to either select a ROI or not. 
        end
    end 
end

