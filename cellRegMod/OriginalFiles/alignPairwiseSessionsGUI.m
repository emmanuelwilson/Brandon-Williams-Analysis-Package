%%Will run CellReg through all session pairs and map common cells
%INPUT: Path to the folder containing subfolders which contain your "ms"
%structures
%
%OUTPUT: 
%   Creates results folder containing all CellReg figures and matrices for
%   each session pair. 
function alignPairwiseSessionsGUI(folder,technique,micperpix, alignment,maxdist,maxrot,rotationsmooth)
%     mice = dir(folder);                                                    %List folder contents
%     mice = {mice(3:end).name};                                             
%     for m = mice
%         mp = [folder '\' m{1}];                                            
        sessions = dir(folder);
        sessions = {sessions(4:end).name};
        
%         % order sessions by date
%         clear dates
        for j = 1:length(sessions)
            seshnums(j) = str2num(sessions{j}(3:end-4));
        end
        [a b] = sort(seshnums);
        sessions = sessions(b);
                
        alignmentMap = repmat({[]},[length(sessions) length(sessions)]);
        for j = 1:length(sessions)
            for k = j+1:length(sessions)
            
                if j+1 > length(sessions) %%% In case odd number of sessions
                    continue
                end
                
                for si = [j k]
                    ref = load([folder '/' sessions{si}]);                
                    [prepped,cellmap] = msExtractSFPsCellReg(ref.ms);
                    if ~isempty(find(sum(sum(prepped,3),2)==0))
                        beep
                        pause
                    end
                    outP = [folder,'/','SegmentsForAlignment/' sessions{si}(1:end-4)];
                    outP2 = [folder,'\Results\',num2str(j),'_',num2str(k)];
                    checkP(outP);
                    checkP(outP2);
                    save(outP,'prepped');
                    if ~isempty(find(diff(cellmap) > 1 ))
                        save(outP2,'cellmap');
                    end
                end

                map = registerCellsGUI([folder, '\','SegmentsForAlignment'], [num2str(j),'_',num2str(k)],technique,micperpix, alignment,maxdist, maxrot,rotationsmooth);
                close all
                copyfile([folder, '\','SegmentsForAlignment\Results\', num2str(j),'_',num2str(k)],outP2);
                rmdir([folder,'\SegmentsForAlignment'],'s');                
                drawnow

%                 isReg = map(all(map~=0,2),:);
                alignmentMap{j,k} = map;%isReg;
            end
        end
        Singlemap = ReorganizeAlignmentMap(alignmentMap);
        save([folder,'\','aligmentMaps'],'alignmentMap','Singlemap','-v7.3');
        
%     end
end