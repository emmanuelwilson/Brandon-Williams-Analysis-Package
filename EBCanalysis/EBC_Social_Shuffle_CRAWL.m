%% Crawl script for EBC shuffling
%Must have all wanted folders under your "Current Folder" 
clear
oldcd = pwd;
matfiles = dir('*.mat');
dirfiles = dir(pwd);
for i = 1 : length(matfiles)    
    for j = 1 : length(dirfiles)
        if dirfiles(j).isdir && contains(dirfiles(j).name,matfiles(i).name(1:end-4))
            load(matfiles(i).name)
            matfiles(i).name
            cd([pwd '/' dirfiles(j).name])
            if exist('normNeuC') && length(dir(pwd)) > 2 
                try
                    load([matfiles(i).name(1:end-4) 'Shuffle.mat'])
                end
                try
                    load([matfiles(i).name(1:end-4) 'Suffle.mat'])
                end
                if ~exist('out')
                    EBC_Social_Shuffle(normNeuC,Angle,DistanceH,1,3,8,[matfiles(i).name(1:end-4) 'Shuffle']);
                end                
            end
            clearvars -except i j matfiles dirfiles oldcd
            cd(oldcd)
            break
        end
    end
end