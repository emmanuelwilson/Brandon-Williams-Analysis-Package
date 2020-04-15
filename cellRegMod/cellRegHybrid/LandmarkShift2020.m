%Find landmark shift between adjacent sessions and saves SFP outlines
%
%INPUT: 
%   -folderpath: folder path containing sessions that are to be compared
%OUPUT:
%   -Shift: Matrix containing FOV shift based off of landmarks. First
%   column represents the horizontal shift, second column represents
%   vertical shift.
%   -Wil overwrite miniscope structure 

function Shift = LandmarkShift2020(folderpath)

folder = dir(folderpath);                                                   %List folder contents
sessions = {folder(3:end).name};                                            %Seperate folder contents
sessions = sessions(3:end);
Shift = [];
count = 0;
msflag = 0;
calflag = 0;
tempsessions = sessions;
for i = 1 : length(sessions)
    if contains(sessions{i},'msDeconvolved')
        if str2num(sessions{i}(14:end-3)) ~= i
            tempsessions{str2num(sessions{i}(14:end-3))} = sessions{i};
        end
    elseif contains(sessions{i},'ms')
        if str2num(sessions{i}(3:end-3)) ~= i
            tempsessions{str2num(sessions{i}(3:end-3))} = sessions{i};
        end
    end
end
sessions = tempsessions;
for i = 1 : length(sessions)-1
    if isempty(sessions(i))                                                %if contents empty skip itteration
        continue
    else
        count = count + 1;
    end
    
    ref = load([folderpath '/' sessions{i}]);                               %Load reference session
    move = load([folderpath '/' sessions{i+1}]);                            %load comparision session
    
    if isfield(ref, 'ms') 
        ref = ref.ms;
        if ~isfield(ref,'outlines')                                        %generate and save footprint outlines 
            ref = SPFoutline(ref);
            ms = ref;
            save([folderpath,'/',sessions{i}],'ms')
        end
    elseif isfield(ref, 'calcium')
        ref = ref.calcium;
        if ~isfield(ref,'outlines')                                        %generate and save footprint outlines
            ref = SPFoutline(ref);
            calcium = ref;
            save([folderpath,'/',sessions{i}],'calcium')
        end
    else
        continue
    end
    if isfield(move,'ms')
        move = move.ms;
        if ~isfield(move,'outlines')
            move = SPFoutline(move);
            ms = move;
            save([folderpath,'/',sessions{i+1}],'ms')
            msflag = 1;
        end
    elseif isfield(move,'calcium')
        move = move.calcium;
        if ~isfield(move,'outlines')
            move = SPFoutline(move);
            calcium = move;
            save([folderpath,'/',sessions{i+1}],'calcium')
        end
    end
    %Find FOV shift between unmodified FOV's
    [wShift, hShift] = msAlignBetweenSessions2019(ref, move,Shift);              %Align sessions
    Shift(count,1) = wShift;
    Shift(count,2) = hShift;
    %Apply previous session shifts
    Xshift(count) = Shift(count,1);
    Yshift(count) = Shift(count,2);
    if count > 1
        Xshift(count) = Xshift(count-1) + Shift(count,1);
        Yshift(count) = Yshift(count-1) + Shift(count,2);
    end
    s = [Xshift(count), Yshift(count)];
    %Modify Footprints with shift
    if msflag
        ms = SFPshift(move,s);
        save([folderpath,'/',sessions{i+1}],'ms')
    elseif calflag
        calcium = SFPshift(move,s);
        save([folderpath,'/',sessions{i+1}],'calcium')
    end
        
end

end