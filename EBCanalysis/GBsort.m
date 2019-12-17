function [] = GBsort(passed, partpassed)
countg = 1;
countb = 1;
goodstop = 1;
badstop = 1;
searchfolder = uigetdir(pwd,'Search folder');
savefolder = uigetdir(pwd,'Save folder');

olddir = pwd;
cd(searchfolder)
load('EBCstats')
folder = dir(pwd);
names = {folder.name};

for i = 3 : length(names)
    if goodstop && ~isempty(find(passed == str2double(names{i}(1:end-7)),1))
        ind = find(passed == str2double(names{i}(1:end-7)),1);
        copyfile([names{i}(1:end-7),'EBC.jpg'],savefolder);        
        movefile([savefolder,'/',names{i}(1:end-7),'EBC.jpg'],[savefolder,'/GOOD_',num2str(out.mrall(passed(ind))),'_Cell',num2str(passed(ind)),'.jpg']);
        countg = countg + 1;
    elseif badstop && ~isempty(find(partpassed == str2double(names{i}(1:end-7)),1))
        ind = find(partpassed == str2double(names{i}(1:end-7)),1);
        copyfile([names{i}(1:end-7),'EBC.jpg'],savefolder);
        movefile([savefolder,'/',names{i}(1:end-7),'EBC.jpg'],[savefolder,'/BAD_',num2str(out.mrall(partpassed(ind))),'_Cell',num2str(partpassed(ind)),'.jpg']);
        countb = countb + 1;
    end
    if (countg > length(passed) || countg > 6) && (countb > length(partpassed) || countb > 6)
        break;
    end
    if badstop && (countb > length(partpassed) || countb > 6)
        badstop = 0;
    end
    if goodstop && (countg > length(passed) || countg > 6)
        goodstop = 0;
    end
end
cd(olddir)
end