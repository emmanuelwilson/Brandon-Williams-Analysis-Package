function []  = SocialProximityMaster_Crawl(p)

paths = genpath(p);
if ispc
    folders = strsplit(paths,';')';
else
    folders = strsplit(paths,':')';
end

for i = 1 : length(folders)    
    try
        cd(folders{i})
        SocialProximityFiring_V3_execute
    end
end
end