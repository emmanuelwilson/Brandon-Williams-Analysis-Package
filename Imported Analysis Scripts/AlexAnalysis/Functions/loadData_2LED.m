function resultingP = loadData_2LED(inFolder,outFolder,envSize,doWindowed)

    if nargin < 4 || isempty(doWindowed);
        doWindowed = false;
    end
        

    warning off all
    if isempty(gcp)
        parpool('local',7);
    end
    pctRunOnAll warning off all

    clc

    allPaths = getFilePaths(inFolder,'timestamp.dat');
%     allPaths = allPaths(end);
%     allPaths = flipud(allPaths);
    resultingP = [];
    for i = 1:length(allPaths)       
        slashInds = find(ismember(allPaths{i},'/'));
        tmp = allPaths{i}(1:find(ismember(allPaths{i},'/'),1,'last')-1);
        fprintf(['\n\t\t\tMouse:  ' tmp '\n']);
        outP = [outFolder '/' tmp(slashInds(1)+1:end) '.mat'];
        
        
        if exist(outP,'file')==2
            fprintf(['\n\t\t\tMouse:  ' tmp(slashInds(2)+1:slashInds(3)-1) '  Already loaded.\n'])
            continue
        end
        
        s = struct;
% 
        fprintf('\n\t*************************************');
        fprintf(['\n\t\t\tMouse:  ' tmp(slashInds(2)+1:slashInds(4)-1)])
        fprintf('\n\t*************************************');
        

        s.properties.name = tmp(slashInds(1)+1:slashInds(2)-1);
        s.properties.session = tmp(slashInds(2)+1:slashInds(3)-1);
        s.properties.trial = tmp(slashInds(3)+1:end);

        fprintf(['\n\n\tTracking position data from:  ' allPaths{i}(1:find(ismember(allPaths{i},'/'),1,'last')-1)])    
        [s.pos.mp s.pos.muninterp] = trackPosition_2LED(allPaths{i}(1:find(ismember(allPaths{i},'/'),1,'last')-1));

        fprintf(['\n\n\tProcessing calcium from:  ' allPaths{i}(1:find(ismember(allPaths{i},'/'),1,'last')-1) '\n'])
        ms = msRunUpdated(allPaths{i}(1:find(ismember(allPaths{i},'/'),1,'last')-1),doWindowed);
        s.calcium = ms;

        close all
        drawnow

        s.environment.type = 'box';
        s.environment.size = envSize;

        s = alignTraceData(allPaths{i}(1:find(ismember(allPaths{i},'/'),1,'last')-1),s);

        outP = [outFolder '/' tmp(slashInds(1)+1:end)];
        checkP(outP);
        save(outP,'-struct','s','-v7.3')
        
        normTrace({outP});
        
        resultingP = [resultingP; {outP}];
    end
end