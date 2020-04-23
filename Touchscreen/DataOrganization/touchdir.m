%% Makes a file directory given the names extrated from the touchscreen software

function [] = touchdir(path)
if isempty(path)
    path = uigetdir();
end
oldcd = pwd;
cd(path);

folder = dir(pwd);
for i = 1 : length(folder)
    if contains(folder(i).name,'.csv')
        fid = fopen('FEB02_100203.csv');
        metanames = textscan(fid,'%s%s%f','delimiter',',');
        fclose(fid);
        
        if length(metanames{1}) == 15 || length(metanames{1}) == 16
                schedulename = metanames{2}{2};
                version = metanames{2}{3};
                enviroment = metanames{2}{4};
                Machinename = metanames{2}{5};
                Datetime = metanames{2}{6};
                Database = metanames{2}{7};
                scheduleRunID = metanames{2}{8};
                Fininalized = metanames{2}{9};
                RecordCount = metanames{2}{10};
                AnimalID = metanames{2}{11};
                GroupID = metanames{2}{12};
                Max_Number_Trials = metanames{2}{13};
                Max_Schedule_Time = metanames{2}{14};
                if length(metanames{1}) == 16
                    User =  metanames{2}{15};
                end
        elseif length(metanames{1}) == 1
            %No meta data
            %Coco's Naming format
            if contains(folder(i).name,'CAM')
                datadash = strfind(database, '_');
                User = 'CAM';
                AnimalID = folder(i).name(4:datadash(1));
                Datetime = folder(i).name(datadash(2)+1:end-4);
                Schedulename = folder(i).name(datadash(2)+1:end-4);
                %Zeeshans Naming format
            elseif contains(folder(i).name,'Zee')
                datadash = strfind(database, '_');
                Machinename = folder(i).name(1:datadash(1));
                Database = folder(i).name(datadash(1):datadash(3));
                folder(i).name(datadash(2)+1:end-4);
                Schedulename = folder(i).name(datadash(3)+1:datadash(4));
                scheduleRunID = folder(i).name(datadash(4):end-4);
                %Andrés Naming format
            else
                datadash = strfind(database, '_');
                AnimalID = folder(i).name(datadash(1):end-4);
                Datetime = folder(i).name(1: datadash(1)-1);
            end

        else 
            warning('Odd number of meta data inputs')
        end

databasepath = database;
databasepath(datadash) = '/';

schedulepath = schedulename;
schedulepath(scheduledash) = '/';

path = [ databasepath, schedulepath, '/', recordingdate];

mkdir(path)
    end
end
    
    %             database = 'User_Region_Task';
%             schedulename = 'TUNL Mouse Exp 1 Stage 1 S3TTL D6sec';
%             recordingdate = 'Date';
%             scheduledash = strfind(schedulename, ' ');
%             datadash = strfind(database, '_');
% if contains(schedulename,'Exp')
%     exploc = strfind(schedulename,'Exp');
%     dashmod = scheduledash;
%     dashmod(dashmod<exploc) = 0;
%     scheduledash(find(dashmod,1)) = [];
% end
% 
% if contains(schedulename,'Stage')
%     exploc = strfind(schedulename,'Stage');
%     dashmod = scheduledash;
%     dashmod(dashmod<exploc) = 0;
%     scheduledash(find(dashmod,1)) = [];
% end
%   