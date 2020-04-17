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
        
        if length(metanames{1}) == 15            
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
        elseif length(metanames{1}) == 1
            database = 'User_Region_Task';
            schedulename = 'TUNL Mouse Exp 1 Stage 1 S3TTL D6sec';
            recordingdate = 'Date';
            scheduledash = strfind(schedulename, ' ');
            datadash = strfind(database, '_');
        else 
            warning('Odd number of meta data inputs')
        end

if contains(schedulename,'Exp')
    exploc = strfind(schedulename,'Exp');
    dashmod = scheduledash;
    dashmod(dashmod<exploc) = 0;
    scheduledash(find(dashmod,1)) = [];
end

if contains(schedulename,'Stage')
    exploc = strfind(schedulename,'Stage');
    dashmod = scheduledash;
    dashmod(dashmod<exploc) = 0;
    scheduledash(find(dashmod,1)) = [];
end
  
databasepath = database;
databasepath(datadash) = '/';

schedulepath = schedulename;
schedulepath(scheduledash) = '/';

path = [ databasepath, schedulepath, '/', recordingdate];

mkdir(path)
end