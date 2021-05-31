function []  = ObjectProximityMaster_Crawl(p)

paths = genpath(p);
if ispc
    folders = strsplit(paths,';')';
else
    folders = strsplit(paths,':')';
end

for i = 1 : length(folders)    
    try
        cd(folders{i})
        ObjectProximity_execute
    end
end
end