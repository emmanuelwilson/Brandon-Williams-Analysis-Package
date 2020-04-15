%% Makes a file directory given the names extrated from the touchscreen software

function [] = touchdir()
folder = dir(pwd);
for i = 1 : length(folder)
    if contains(folder(i).name,'.csv')
        readMatrix([folder(i).folder,'/',folder(i).name])
database = 'User_Region_Task';

schedulename = 'TUNL Mouse Exp 1 Stage 1 S3TTL D6sec';

recordingdate = 'Date';

scheduledash = findstr(schedulename, ' ');

datadash = findstr(database, '_');

if contains(schedulename,'Exp')
    exploc = findstr(schedulename,'Exp');
    dashmod = scheduledash;
    dashmod(dashmod<exploc) = 0;
    scheduledash(find(dashmod,1)) = [];
end

if contains(schedulename,'Stage')
    exploc = findstr(schedulename,'Stage');
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