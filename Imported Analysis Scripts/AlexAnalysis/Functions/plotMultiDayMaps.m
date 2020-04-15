function plotMultiDayMaps(paths)
    
    warning off all
    
    clc
    fprintf('\n')
    %%% Reliability constaint
    
    %% Split by animal
    piece = [];
    ag = [];
    spiece = [];
    for i = 1:length(paths)
        ind = find(ismember(paths{i},'/'),1,'last')-1;
        piece = [piece; {paths{i}(1:ind)}];
        spiece = [spiece; {paths{i}(ind+2:end-4)}];
    end
    upiece = unique(piece);
    
    for mi = 1:length(upiece)
        fprintf(['\n\tMouse:  ' num2str(upiece{mi}) '\n'])
        isM = find(ismember(piece,upiece(mi)));
        sessions = paths(isM);
        s = load(paths{isM(1)});
        alignMap = s.alignment.allAlignmentMap;
        am = [];
        for si = 1:length(sessions)
            s = load(sessions{si});  
            [isIn isMostRecent indexSinceIn] = isInROI(s.processed.p,s.processed.roi.door);
            [isInRoom] = isInROI(s.processed.p,s.processed.roi.room);
            allMasks = [{isMostRecent(1,isInRoom)}  {isMostRecent(2,isInRoom)}];
            gT = s.processed.trace(alignMap(:,si),:);
            mapA = mkTraceMaps(s.processed.p(:,isInRoom),gT(:,isInRoom),allMasks{1},[13 13]);
            mapB = mkTraceMaps(s.processed.p(:,isInRoom),gT(:,isInRoom),allMasks{2},[13 13]);
            map2 = cat(1,mapA,nan([1 length(mapA(1,:,1)) length(mapA(1,1,:))]),mapB);
            am = cat(2,am,nan(length(map2(:,1,1)),1,length(map2(1,1,:))),map2);
        end
        
        for k = 1:length(am(1,1,:))
            tmp = am(:,:,k);
            figure(1)
            set(gcf,'position',[50 50 10.*length(tmp(1,:)) 10.*length(tmp(:,1))])
            imagesc(tmp)
            alpha(double(~isnan(tmp)))
            axis off
            colormap jet
            
            slashInds = find(ismember(sessions{1},'/'));
            outP = ['Plots/AlignedCellMaps/' sessions{1}(slashInds+1:end-4) '_Cell_' num2str(k)];
            saveFig(gcf,outP,[{'tiff'} {'pdf'}])
            close all
            drawnow
        end
        
    end
end