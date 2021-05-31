%% Super session creator: will take all calcium videos given path and sort them in order for analysis
%   INPUT: path containing subfolders "Habituation", "trial1", "trial2" and "trial3"

function [] = SuperSessionSITv2(p)

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
prevfolder = [];
for i = 1 : length(folders)
    if ~isempty(folders{i})
        d = dir(folders{i});
        fnames = {d.name};
        if i == 1 %~isempty(find(strcmp(fnames,'ms.mat'),1))
            subfolders = strsplit(folders{i},'\');
            subfolders = strsplit(folders{i},'\');
            newName = subfolders{1};
            if isempty(find(strcmp(newName, mouse),1))
                mouse = newName;
                mkdir('ConcactenatedSession')
            end
        end
        j = [];
        if ~isempty(find(strcmp(fnames,'ms.mat'),1))
            subfolders = strsplit(folders{i},'\');
            subfolders = strsplit(folders{i},'\');
            j = length(subfolders);
            if ~isempty(find(contains(subfolders{j}, 'trial'),1))                
                    cd(folders{i})
                    mscamdir = dir('msCam*.avi');
                    catvidcount = dir('../../ConcactenatedSession/*.avi');
                    frameMapcount = dir('../../ConcactenatedSession/frameMap*.mat');
                    badframescount = dir('../../ConcactenatedSession/badframes*.mat');
                    HDcount = dir('../../ConcactenatedSession/HeadTrackingData*.mat');
                    for vid = 1 : length(mscamdir)
                        copyfile(['msCam' num2str(vid) '.avi'], ['../../ConcactenatedSession/msCam' num2str(length(catvidcount) + vid) '.avi']);
                    end
                    copyfile('frameMap.mat', ['../../ConcactenatedSession/frameMap' num2str(length(frameMapcount) + 1) '.mat']);
                    copyfile('badframes.mat', ['../../ConcactenatedSession/badframes' num2str(length(badframescount) + 1) '.mat']);
                    copyfile('HeadTrackingData.mat', ['../../ConcactenatedSession/HeadTrackingData' num2str(length(HDcount) + 1) '.mat']);
%                     break                
            elseif ~isempty(find(contains(subfolders{j}, 'habituation'),1))
                cd(folders{i})
                mscamdir = dir('msCam*.avi');
                catvidcount = dir('../../ConcactenatedSession/*.avi');
                frameMapcount = dir('../../ConcactenatedSession/frameMap*.mat');
                badframescount = dir('../../ConcactenatedSession/badframes*.mat');
                HDcount = dir('../../ConcactenatedSession/HeadTrackingData*.mat');
                %                     if ~isempty(catvidcount)
                %                         error('out of sequence')
                %                     end
                for vid = 1 : length(mscamdir)
                    copyfile(['msCam' num2str(vid) '.avi'], ['../../ConcactenatedSession/msCam' num2str(length(catvidcount) + vid) '.avi']);
                end
                copyfile('frameMap.mat', ['../../ConcactenatedSession/frameMap' num2str(length(frameMapcount) + 1) '.mat']);
                copyfile('badframes.mat', ['../../ConcactenatedSession/badframes' num2str(length(badframescount) + 1) '.mat']);
                copyfile('HeadTrackingData.mat', ['../../ConcactenatedSession/HeadTrackingData' num2str(length(HDcount) + 1) '.mat']);
%                 break
            end
            
            
        end
    end
end