%%Will run CellReg through all session pairs and map common cells
%INPUT: Path to the folder containing subfolders which contain your "ms"
%structures
%
%OUTPUT: 
%   Creates results folder containing all CellReg figures and matrices for
%   each session pair. 
function alignPairwiseSessions(folder)                                          
        sessions = dir(folder);
        sessions = {sessions(5:end).name};
        
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
                    if (cellmap(1) > 1) || ~isempty(find(diff(cellmap) > 1 ))
                        save([folder,'\Results\cellmap', num2str(si),'.mat'],'cellmap');
                    end
                end
                map = registerCellsNoGUI([folder, '\','SegmentsForAlignment'], [num2str(j),'_',num2str(k)], j , k);                
                close all
                copyfile([folder, '\','SegmentsForAlignment\Results\', num2str(j),'_',num2str(k)],outP2);
                rmdir([folder,'\SegmentsForAlignment'],'s');                
                drawnow                

%                 isReg = map(all(map~=0,2),:);
                alignmentMap{j,k} = map;%isReg;
            end
        end
        
        cormat = cormatrix([folder,'\Results']);
        save([folder,'\cormat.mat'],'cormat');
        
        Singlemap = ReorganizeAlignmentMap(alignmentMap);
        probmap = RegProbMap([folder,'\Results']);
        [Singlemap, avg_psame] = ElimConflict(Singlemap,alignmentMap,probmap);
        Singlemap = FindMissingCells(Singlemap,folder);
        save([folder,'\','aligmentMaps'],'alignmentMap','Singlemap','probmap','avg_psame','-v7.3');
end