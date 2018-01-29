function frame = msReadFrame(ms,frameNum,columnCorrect, align, dFF)
%this function reads in variables from ms structure and outputs a frame
%from a video object within the structure. 
    
    vidNum = ms.vidNum(frameNum);                           % here we call in a video from a vidObj,
    vidFrameNum = ms.frameNum(frameNum);                    % pick a frame given frameNUm, and from that 
    frame = double(ms.vidObj{vidNum}.read(vidFrameNum));    % video store in variable frame 
    
%-------------------- Use same processing algorithm as in msAutoSegment.m to increase SNR ---------------------------------- 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------------------------------------------------------------
    
    %The next section of code would compare multiple frames together and
    %shift them to compare statistically the error between them. 

    if (align)
        frame = frame(((max(ms.hShift(:,ms.selectedAlignment))+1):(end+min(ms.hShift(:,ms.selectedAlignment))-1))-ms.hShift(frameNum,ms.selectedAlignment), ...
                      ((max(ms.wShift(:,ms.selectedAlignment))+1):(end+min(ms.wShift(:,ms.selectedAlignment))-1))-ms.wShift(frameNum,ms.selectedAlignment));
    end
    
    if (columnCorrect)
        %frame = frame - ms.columnCorrection + ms.columnCorrectionOffset;
    end
    
    if (dFF)
%         idx = ms.minFrame{ms.selectedAlignment}<80;
        frame = frame./ms.minFrame{ms.selectedAlignment}-1;
%         frame = frame./ms.meanFrame{ms.selectedAlignment}-1;
%         frame(idx) = 0;
    end
    
end


%a=(ms.vidObj{vidNum}.read(vidFrameNum))
%mask = roipoly(a)
%imshow(im2double(mask).*im2double(a))
%s=sum(mask(:))