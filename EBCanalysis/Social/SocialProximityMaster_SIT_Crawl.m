function []  = SocialProximityMaster_SIT_Crawl(p)

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
        if ~isempty(find(strncmp(fnames,'frameMap.mat',10),1)) && ~isempty(find(strncmp(fnames,'SITstartFrame.mat',10),1))
            try
                cd(folders{i})
                SocialProximityFiring_SIT_execute
            end
        end
    end
end
end