function batchMerge(inFolder,outFolder)

    tf = getFilePaths(inFolder,'timestamp.dat');
    piece = [];
    snames = [];
    for i = 1:length(tf)
        sinds = find(ismember(tf{i},'/'));
        t = tf{i}(sinds(3)+1:sinds(4)-1);
        piece = [piece; {[tf{i}(sinds(2)+1:sinds(3)-1) '/' t(1:find(ismember(t,'_'),1,'last')-1)]}];
    end
    
    clc
    fprintf('\n\tMerging Sessions for:\n')
    upiece = unique(piece);
    for mi = 1:length(upiece)
        fprintf(['\n\t\t' upiece{mi}])
        isS = ismember(piece,upiece(mi));
%         if sum(isS)~=2
%             continue
%         end
        catMiniscopeFiles(tf(isS),[outFolder '/' upiece{mi}]);
    end    
end