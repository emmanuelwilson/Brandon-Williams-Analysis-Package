%% Will convert # of frames to seconds and output cell array for graphing label
%INPUT:
%   -L: length of array in frames 
%   -framerate: how many frames in a second
%OUTPUT:
%   -timeLabels: cell array containing characters of wanted time labels in
%   seconds.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Emmanuel Wilson

function timeLabels = Frame2SecLabels(L,framerate,Ticrate)

Ttot = floor(L/framerate);
for i = 0 : Ticrate: Ttot
    timeLabels{i+1} = num2str(i);
end

end