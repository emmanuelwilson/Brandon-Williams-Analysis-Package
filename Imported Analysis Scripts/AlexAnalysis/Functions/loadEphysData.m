function loadEphysData(inFolder,outFolder,overwrite)

    addpath(genpath('Z:\MarkAnalysis'));

    allPaths = getFilePaths(inFolder,'_object.mat');
    for i = 1:length(allPaths)       
        slashInds = find(ismember(allPaths{i},'/'));
        tmp = allPaths{i}(1:find(ismember(allPaths{i},'/'),1,'last')-1);
        fprintf(['\n\t\t\tMouse:  ' tmp '\n']);
        outP = [outFolder '/' tmp(slashInds(1)+1:end) '.mat'];
        
        
        if exist(outP,'file')==2
            fprintf(['\n\t\t\tMouse:  ' tmp(slashInds(2)+1:slashInds(3)-1) '  Already loaded.\n'])
            continue
        end
%       
%         if ~ismember({tmp(slashInds(2)+1:slashInds(3)-1)},[{'AKCA105'}])
%             continue
%         end
        
        s = struct;
% 
        fprintf('\n\t*************************************');
        fprintf(['\n\t\t\tMouse:  ' tmp(slashInds(2)+1:slashInds(4)-1)])
        fprintf('\n\t*************************************');
        

        fprintf(['\n\n\t\tLoading Ephys:  ' tmp(slashInds(2)+1:slashInds(4)-1) ' (' num2str(i) ' of ' num2str(length(allPaths))  ')'])

        load([tmp '/_object.mat']);
        
        s.properties.rat = tmp(slashInds(end-2)+1:slashInds(end-1)-1);
        s.properties.date = tmp(slashInds(end-1)+1:slashInds(end-1)+8);


        s.pos.p = [root.x'; root.y'] .* 0.21;
        s.pos.ts = root.ts';
        s.pos.hd = deg2rad(root.headdir');
        s.pos.vel = root.vel';
        s.environment.size = [36 36]; % [66 42]
        
        s.p.cm2pix = 2.5;
        s.p.kern_sd = 4;

        s.unit = repmat({[]},size(root.spike));
        s.processed.trace = zeros(nansum(~cellfun(@isempty,s.unit(:))),length(s.pos.p));
        count = 0;
        for t = 1:length(s.unit(:,1))
            for c = 1:length(s.unit(1,:))
                if ~isempty(root.spike(t,c).ts)
                    count = count+1;
                    s.unit{t,c} = root.spike(t,c).ts;
                    
                    ts = allign(s.pos.ts,root.spike(t,c).ts');
                    s.processed.trace(count,ts) = 1;
                end
            end
        end
        s.processed.p = s.pos.p;
        
        checkP(outP);
        save(outP,'-struct','s','-v7.3')
    end
    rmpath(genpath('Z:\MarkAnalysis'));
end