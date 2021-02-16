%% Execute automated CellReg2020: Pairwise Comparison

%Select the folder which contains all trail information
folderpath = uigetdir(pwd,'Select Folder Containing all sessions miniscope variables');
oldcd = pwd;
cd(folderpath)
dir(pwd);
% mkdir(folderpath, 'OriginalFiles')
% copyfile(folderpath,[folderpath,'/OriginalFiles']);
folder = dir(folderpath);

%Find FOV shift between neighboring sessions
[wShift,hShift,sessions,nonRigid,combs] = Automation_LandmarkShift2020_Pairwise(folderpath);

%Align sessions using CellReg
alignNwiseSessions2020_Auto_Pairwise(wShift,hShift,sessions,combs,folderpath,nonRigid)
