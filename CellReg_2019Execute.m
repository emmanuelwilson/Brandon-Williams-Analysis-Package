%% Execute CellReg2019

%Select the folder which contains all trail information
folderpath = uigetdir(pwd,'Select Folder Containing all sessions miniscope variables');
mkdir(folderpath, 'OriginalFiles')
copyfile(folderpath,[folderpath,'/OriginalFiles']);
folder = dir(folderpath);

%Find FOV shift between neighboring sessions
Shift = LandmarkShift(folderpath);

%Align sessions using CellReg
alignPairwiseSessions(folderpath)

cormat = cormatrix([folderpath,'\Results']);
save([folderpath,'\cormat.mat'],'cormat');
