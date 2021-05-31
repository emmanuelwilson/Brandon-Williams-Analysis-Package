
function [] = IdentifyBadFrames_crawl(p)

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
        if ~isempty(find(strncmp(fnames,'HeadTrackingData.mat',16),1)) && isempty(find(strncmp(fnames,'badframes.mat',10),1))
            cd(folders{i});                                  %Change current folder            
                load('HeadTrackingData.mat')
                figure(1)
                plot(SINKdata(:,1))                
                hold on
                plot(SINKdata(:,2))
                
                t = input('Up to which frame is bad?');
                save('badframes','t');
                
                close all
        end
    end
end
