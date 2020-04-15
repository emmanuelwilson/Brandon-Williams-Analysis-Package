function selectPoints(paths,overwrite)
    
    clc
    fprintf('\nSelect port locations of interest:\n')
    
    if nargin < 2
        overwrite = false;
    end

    for i = 1:length(paths)
        
        fprintf(['\n\t' paths{i}]);
        
        s = load(paths{i});
        
        if ~overwrite && isfield(s.processed,'locs')
            continue
        end
        
%         figure(1)
%         set(gcf,'position',[50 50 800 800])
%         plot(s.processed.p(1,:),s.processed.p(2,:));
%         axis square
%         axis equal
%         [x y] = getpts();
%         s.processed.locs = [x'; y']';
%         close all
%         drawnow;
        
        es = s.processed.envSize;
        s.processed.locs = [es./2; es(1) 0; 0 es(2); 0 0; es];
    
        save(paths{i},'-struct','s','-v7.3');
    end
end