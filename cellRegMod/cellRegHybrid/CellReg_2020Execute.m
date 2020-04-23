%% Execute CellReg2020

%Select the folder which contains all trail information
folderpath = uigetdir(pwd,'Select Folder Containing all sessions miniscope variables');
oldcd = pwd;
cd(folderpath)
dir(pwd);
mkdir(folderpath, 'OriginalFiles')
copyfile(folderpath,[folderpath,'/OriginalFiles']);
folder = dir(folderpath);

%Find FOV shift between neighboring sessions
Shift = LandmarkShift2020(folderpath);

%Align sessions using CellReg
alignNwiseSessions2020(folderpath)
