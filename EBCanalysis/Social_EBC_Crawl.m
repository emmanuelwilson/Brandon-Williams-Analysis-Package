%% Run Social EBC on all Dataset variables in folder

function [] = Social_EBC_Crawl(p)
oldcd = pwd;
cd(p)
matfiles = dir('*.mat');
dirfiles = dir(pwd);

for i = 6 : length(matfiles)
    predone = 0;
    for j = 1 : length(dirfiles)
        if dirfiles(j).isdir && contains(dirfiles(j).name,matfiles(i).name(1:end-4))
            predone = 1;
            break
        end
    end
    if ~predone
        load(matfiles(i).name)
        matfiles(i).name
        if exist('normNeuC')
            msEgoCentricRateMapSplitEvenOddParallel_Social(normNeuC,Angle,DistanceH,0,[1:length(NeuC(:,1))],3,8,matfiles(i).name(1:end-4));
            close all
            cd(p)
        end
    end
end
cd(oldcd);
end