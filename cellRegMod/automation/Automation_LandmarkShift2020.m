%Find landmark shift between adjacent sessions and saves SFP outlines
%
%INPUT: 
%   -folderpath: folder path containing sessions that are to be compared
%OUPUT:
%   -Shift: Matrix containing FOV shift based off of landmarks. First
%   column represents the horizontal shift, second column represents
%   vertical shift.
%   -Wil overwrite miniscope structure 

function [wShift,hShift,sessions,shifts,nonRigid,combs] = Automation_LandmarkShift2020(folderpath)

prompt = 'Non-Rigid Registration?:Y/N ';
str = input(prompt,'s');
if isempty(str) || str == 'Y' || str == 'y'
    nonRigid = true;
else
    nonRigid = false;
end

folder = dir(folderpath);                                                   %List folder contents
sessions = {folder.name};                                            %Seperate folder contents
nFold = 2;
%make sure that the sessions to be analyzed are not in the these folders
for i = length(folder):-1:1
    if folder(i).isdir
        sessions(i) = [];
    end
end

Shift = [];
count = 0;
msflag = 0;
calflag = 0;
%sort the sessions in numbered order
tempsessions = sessions;
for i = 1 : length(sessions)
    if ~isempty(str2num(sessions{i}(end-5:end-3)))        
        tempsessions{str2num(sessions{i}(end-5:end-3))} = sessions{i};
    elseif ~isempty(str2num(sessions{i}(end-4:end-3)))        
        tempsessions{str2num(sessions{i}(end-4:end-3))} = sessions{i};
    end
end
%save all of the sessions in one variable location
sessions = cell(length(tempsessions),1);
for i = 1 : length(tempsessions)
    try
        sessions{i} = load([folderpath '/' tempsessions{i}],'calcium');
        try            
            sessions{i}.processed = load([folderpath '/' tempsessions{i}],'processed');
            if isempty(fieldnames(sessions{i}.processed))
                sessions{i} = rmfield(sessions{i},'processed');
            end
        end
    catch
        sessions{i} = load([folderpath '/' tempsessions{i}],'ms');
        calcium = sessions{i};
        sessions{i} = calcium;
    end
    if isempty(fieldnames(sessions{i}))
        sessions{i} = load([folderpath '/' tempsessions{i}],'ms');
        calcium.calcium = sessions{i}.ms;
        sessions{i} = calcium;
    end
    try
        sessions{i}.calcium  = rmfield(sessions{i}.calcium ,'numFiles');
        sessions{i}.calcium  = rmfield(sessions{i}.calcium ,'numFrames');
        sessions{i}.calcium  = rmfield(sessions{i}.calcium ,'vidNum');
        sessions{i}.calcium  = rmfield(sessions{i}.calcium ,'vidObj');
        sessions{i}.calcium  = rmfield(sessions{i}.calcium ,'frameNum');
        sessions{i}.calcium  = rmfield(sessions{i}.calcium ,'maxFramesPerFile');
        sessions{i}.calcium  = rmfield(sessions{i}.calcium ,'dateNum');
        sessions{i}.calcium  = rmfield(sessions{i}.calcium ,'analysis_time');
        sessions{i}.calcium  = rmfield(sessions{i}.calcium ,'ds');
        sessions{i}.calcium  = rmfield(sessions{i}.calcium ,'shifts');
        sessions{i}.calcium  = rmfield(sessions{i}.calcium ,'Options');
        sessions{i}.calcium  = rmfield(sessions{i}.calcium ,'Centroids');
        sessions{i}.calcium  = rmfield(sessions{i}.calcium ,'CorrProj');
        sessions{i}.calcium  = rmfield(sessions{i}.calcium ,'PeakToNoiseProj');
        sessions{i}.calcium  = rmfield(sessions{i}.calcium ,'FiltTraces');
        sessions{i}.calcium  = rmfield(sessions{i}.calcium ,'RawTraces');
        sessions{i}.calcium  = rmfield(sessions{i}.calcium ,'analysis_duration');
    end
        
    sessions{i} = SPFoutline(sessions{i}.calcium);
    sessions{i}.meanFiltFrame = fourrierfilter(sessions{i}.meanFrame,10,70);     
end
combs = nchoosek(1:length(sessions),nFold);

[wShift,hShift] = msAlignBetweenSessions2020(sessions,combs);
shifts = cell(length(sessions),length(sessions));
for i = 1 : length(combs)
    s = sessions{combs(i,2)};
    ms = SFPshift(s,[wShift(combs(i,1),combs(i,2)),hShift(combs(i,1),combs(i,2))]);    
    shifts{combs(i,1),combs(i,2)} = ms.SFPs;
end
save([folderpath,'/'],'shifts')

end