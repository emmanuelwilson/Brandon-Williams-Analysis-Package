function out = loadSocialProximitydata_Amy(p)

paths = genpath(p);
if ispc
    folders = strsplit(paths,';')';
else
    folders = strsplit(paths,':')';
end
socialProxSessions = [];
exclude = [];
mouse = [];
countM = 0;
fristtime = 0;
micenames = [];
for i = 1 : length(folders)
    if ~isempty(folders{i})
        d = dir(folders{i});
        fnames = {d.name};
        if ~isempty(find(strcmp(fnames,'SocialProximity_V2.mat'),1))
            cd(folders{i});                                  %Change current folder
            load('SocialProximity_V2.mat')
            load('ms.mat')
            subfolders = strsplit(folders{i},'\');
            for j = length(subfolders) : -1 : 1
                if ~isempty(find(contains(subfolders{j}, 'camkii-'),1)) && isempty(find(contains(subfolders{j}, 'trial'),1)) && isempty(find(contains(subfolders{j}, 'habituation'),1))
                    newName = subfolders{j};                   
                    if isempty(find(strcmp(newName, mouse),1))
                        countM = countM +1;
                        mouse = newName;
                        micenames{countM} = mouse;
                        socialProxSessions = cat(2,socialProxSessions,cell(4,1));
                        exclude = cat(2,exclude,cell(4,1));
                    end               
                    socialProxSessions{count,countM} = SocialProximity;
                    exclude{count,countM} = ms.exclude;
                    break
                elseif ~isempty(find(contains(subfolders{j}, 'trial'),1))
                    count = str2num(subfolders{j}(end))+1;
                elseif ~isempty(find(contains(subfolders{j}, 'habituation'),1))
                    count = 1;
                end
            end
        end
    end
    out.socialProxSessions = socialProxSessions;
    out.micenames = micenames;
    out.exclude = exclude;
end