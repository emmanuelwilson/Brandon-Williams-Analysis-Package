%Find landmark shift between adjacent sessions and saves SFP outlines
%
%INPUT: 
%   -folderpath: folder path containing sessions that are to be compared
%OUPUT:
%   -Shift: Matrix containing FOV shift based off of landmarks. First
%   column represents the horizontal shift, second column represents
%   vertical shift.
%   -Wil overwrite miniscope structure 
function Shift = LandmarkShift(folderpath)

folder = dir(folderpath);                                                   %List folder contents
sessions = {folder(3:end).name};                                            %Seperate folder contents
sessions = sessions(2:end);
Shift = [];
count = 0;

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
    else
        continue
    end
    if isfield(move,'ms')
        move = move.ms;
        if ~isfield(move,'outlines')
            move = SPFoutline(move);
            ms = move;
            save([folderpath,'/',sessions{i+1}],'ms')
        end
    end
    [wShift, hShift] = msAlignBetweenSessions2018(ref, move,Shift);              %Align sessions        
    Shift(count,1) = wShift;
    Shift(count,2) = hShift;    
    ms = SFPshift(move,Shift(count,:));
    save([folderpath,'/',sessions{i+1}],'ms')
end


for i  = 1 : length(Shift)
    if i == length(Shift)
        continue
    else        
        for j = i+1 : length(Shift)
            if j == i+1 
                Yshift(i,j) = Shift(i,1);
                Xshift(i,j) = Shift(i,2);
            else
                Yshift(i,j) = Yshift(i,j-1) + Shift(i,1);
                Xshift(i,j) = Xshift(i,j-1) + Shift(i,2);
            end
        end
    end
end


end