function out = loadOjbectProximitydata_NOR_Amy(p)

paths = genpath(p);
if ispc
    folders = strsplit(paths,';')';
else
    folders = strsplit(paths,':')';
end
objectProxSessions = [];
exclude = [];
folderpaths = [];
mouse = [];
sub = [];
countM = 0;
fristtime = 0;
micenames = [];
for i = 1 : length(folders)
    if ~isempty(folders{i})
        d = dir(folders{i});
        fnames = {d.name};
        if ~isempty(find(strcmp(fnames,'SocialProximity_V1_2.mat'),1))
            cd(folders{i});                                  %Change current folder
            load('SocialProximity_V1_2.mat')
            load('ms.mat')
            load('ObjectStatsHab.mat')
            subfolders = strsplit(folders{i},'\');
            for j = length(subfolders) : -1 : 1
                if ~isempty(find(contains(subfolders{j}, 'camkii-'),1)) && isempty(find(contains(subfolders{j}, 'trial'),1))
                    newName = subfolders{j};
                    newsub = subfolders{j-1};
                    if isempty(find(strcmp(newName, mouse),1)) || strcmp(newName, mouse) && ~strcmp(newsub, sub)
                        countM = countM +1;
                        mouse = newName;
                        sub = newsub;
                        micenames{countM} = mouse;                        
                        objectProxSessions = cat(2,objectProxSessions,cell(2,1));                        
                        if isempty(exclude)
                            exclude = cell(2,1);
                        else
                            exclude = cat(2,exclude,cell(2,1));
                        end
                        folderpaths= cat(2,folderpaths,cell(2,1));
                    end               
                    objectProxSessions{count,countM} = SocialProximity;
                    if isfield(ms,'exclude')
                        exclude{count,countM} = ms.exclude;
                    end
                    folderpaths{count,countM} = folders{i};
                    break
                elseif ~isempty(find(contains(subfolders{j}, 'trial'),1))
                    count = str2num(subfolders{j}(end));
                end
            end
        end
    end
    out.socialProxSessions = objectProxSessions;
    out.micenames = micenames;
    out.exclude = exclude;
    out.folderpaths= folderpaths;
end