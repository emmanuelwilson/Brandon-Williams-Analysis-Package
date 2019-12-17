function [] = Gsort(passed,searchfolder,savefolder)
if isempty(searchfolder)
    searchfolder = uigetdir(pwd,'Search folder');
end
if isempty(savefolder)
    savefolder = uigetdir(pwd,'Save folder');
end

olddir = pwd;
cd(searchfolder)
load('EBCstats')
folder = dir(pwd);
names = {folder.name};

for i = 3 : length(names)-1
    ind = find(passed == str2double(names{i}(1:end-7)),1);
    if ~isempty(ind)
        copyfile([names{i}(1:end-7),'EBC.jpg'],savefolder);
%     movefile([savefolder,'/',names{i}(1:end-7),'EBC.jpg'],[savefolder,'/GOOD_',num2str(out.mrall(passed(ind))),'_Cell',num2str(passed(ind)),'.jpg']);    
        movefile([savefolder,'/',names{i}(1:end-7),'EBC.jpg'],[savefolder,'/GOOD_Cell',num2str(passed(ind)),'Wall.jpg']);
    end
end
cd(olddir)
end