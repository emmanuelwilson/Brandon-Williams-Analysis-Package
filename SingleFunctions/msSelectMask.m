function vidObj = msSelectMask(vidObj,downSamp)
%MSREMOVECOLUMNVARIATION Generates the vector used to correct for ADC noise
%   downSamp will cut down on the number of frames used in the calculation

    %meanFrame = zeros(vidObj.height,vidObj.width); %allocate memory
    meanFrame = zeros(vidObj.alignedHeight(vidObj.selectedAlignment),vidObj.alignedWidth(vidObj.selectedAlignment)); %allocate memory
    count = 0;
    ask = 0;
    
    for frameNum=1:downSamp:vidObj.numFrames
        count = count + 1;
        
           %----------------------- Zaki ----------------------------------
           if ask ==0
               userInput = upper(input('Define ROI? (Y/N)','s'));
               if(strcmp(userInput,'Y'))
                   selection_frame = msReadFrame(vidObj,frameNum,false,true,false);
                   figure;
                   vidObj.mask = roipoly(uint8(selection_frame));
               else
                   vidObj.mask = ones(size(meanFrame));
               end
                ask = 1;
           end
           %---------------------------------------------------------------
        
        %meanFrame = meanFrame + double(msReadFrame(vidObj,frameNum,false,false,false));
%         meanFrame = meanFrame + double(vidObj.mask).*double(msReadFrame(vidObj,frameNum,false,true,false));
%         if (mod(frameNum,1+100*downSamp)==0)
%             display(['Calculating column correction. ' num2str(frameNum/vidObj.numFrames*100) '% done'])
%         end
%     end
%     
%     % creates correction frame used to remove ADC noise
%     %vidObj.columnCorrection = round(repmat(mean(meanFrame/count,1),vidObj.height,1));
%     vidObj.columnCorrection = round(repmat(mean(meanFrame/count,1),vidObj.alignedHeight(vidObj.selectedAlignment),1));
%     vidObj.columnCorrectionOffset = mean(vidObj.columnCorrection(:));
end
