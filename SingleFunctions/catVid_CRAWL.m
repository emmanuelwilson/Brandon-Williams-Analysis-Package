%% CatVids crawling, will concactenate all videos with a given prefix in each folder within the path

function [] = catVid_CRAWL(p)

paths = genpath(p);
if ispc
    folders = strsplit(paths,';')';
else
    folders = strsplit(paths,':')';
end
for i = 1 : length(folders)
%     folders{i}
    if ~isempty(folders{i})        
        d = dir(folders{i});
        fnames = {d.name};
        if isempty(find(strcmp(fnames,'behavCamCat.avi'),1)) && ~isempty(find(strcmp(fnames,'behavCam1.avi'),1))
            cd(folders{i});                                  %Change current folder
            try
                catVids(pwd, 'behavCam');
            end
        end
    end
end
end