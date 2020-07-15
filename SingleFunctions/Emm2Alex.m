%% Emm2Alex will convert all calcium pipeline "ms" files into Alex's "calcium" format with post-processing and SPFs exclusion
%INPUT: -endPath: Will save results in designated folder with the original
%           file folder location as its new file name. 
%       -startPath: location of folder with subfolders containing the
%           wanted files to convert


function [] = Emm2Alex(endPath,startPath)

paths = genpath(startPath);
folders = strsplit(paths,';')';

for i = 1 : length(folders)
    d = dir(folders{i});
    fnames = {d.name};
    if ~isempty(find(strcmp(fnames,'ms.mat'),1)) && ~isempty(find(strncmp(fnames,'timestamp',1),1))
        load([folders{i} '/ms.mat']);
        s.calcium = ms;
        s = alignTraceData_noTracking(folders{i},s);
        s = Emm_converter_excludeSFPs(folders(i),s);
        if ispc
            varnames = strsplit(folders{i},'\');
        else
            varnames = strsplit(folders{i},'/');
        end
        save([endPath,'/',varnames{end}],'-struct','s','-v7.3');
    end
end