function [] = SocialEBCPrelimenaryFigures_reorganize(p)

paths = genpath(p);
if ispc
    folders = strsplit(paths,';')';
else
    folders = strsplit(paths,':')';
end
for i = 1 : length(folders)
    if ~isempty(folders{i})
        d = dir(folders{i});
        fnames = {d.name};
        if ~isempty(find(strncmp(fnames,'PassedL',1),1)) && ~isempty(find(strncmp(fnames,'PassedR.mat',1),1))
            cd(folders{i});
            subfolders = strsplit(folders{i},'\');
            try
                for j = length(subfolders) : -1 : 1
                    if ~isempty(find(contains(subfolders{j}, 'trial'),1))
                        copyfile('PassedLeft',['E:\SocialPrelimFigures\Trial' subfolders{j}(end) '/Left']);
                        copyfile('PassedRight',['E:\SocialPrelimFigures\Trial' subfolders{j}(end) '/Right']);
                    elseif ~isempty(find(contains(subfolders{j}, 'habituation'),1))
                        copyfile('PassedLeft','E:\SocialPrelimFigures\Habituation/Left');
                        copyfile('PassedRight','E:\SocialPrelimFigures\Habituation/Right');
                    end
                end
            end
        end
    end
end
end