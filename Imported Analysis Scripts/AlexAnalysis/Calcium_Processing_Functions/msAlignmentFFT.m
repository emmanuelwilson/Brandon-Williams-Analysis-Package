function ms = msAlignmentFFT(ms,plotting)
%Takes use of function sbxalgin (a recursive alignment algorithm) to find
%the shift between image frames. This then saves the calculated data into 
%ms.alignmentHeight and ms.alignmentWidth  
    
    for ROINum = 1:size(ms.alignmentROI,2) %runs from 1 to the number of rows in alignmentROI 
%         display(['Alignment ' num2str(ROINum) '/' num2str(size(ms.alignmentROI,2))]); %reads out the alignment 
        rect = ms.alignmentROI(:,ROINum); %saving the first row of alignmentROI in rect.  
        if(mod(round(min(rect([3 4]))),2)==0)
            rect([3 4]) = rect([3 4]) -1;
        end
        ROI = uint16([rect(1) rect(1)+rect(3) rect(2) rect(2)+rect(4)]);
        tic
        fprintf(['\tCalculating shift between frames. Percent frames aligned:      ']);
        r = sbxalign(ms,1:ms.numFrames,ROI,plotting);%uses an alignmnet function 
                                                     %to recursively align
                                                     %images and produce
                                                     %statistics, for which
                                                     %are then saved into
                                                     %variable r. 
        fprintf(['\t']);
        toc
        
        ms.hShift(:,ROINum) = r.T(:,1);%assigns the statistics from sbxalign to hShift 
        ms.wShift(:,ROINum) = r.T(:,2);%assigns the statistics from sbxalign to wShift
        ms.alignedHeight(ROINum) = ms.height - (max(ms.hShift(:,ROINum))- min(ms.hShift(:,ROINum))+1); %
        ms.alignedWidth(ROINum) = ms.width - (max(ms.wShift(:,ROINum))-min(ms.wShift(:,ROINum))+1); %
    end
end

