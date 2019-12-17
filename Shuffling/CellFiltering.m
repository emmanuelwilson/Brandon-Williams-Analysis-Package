function [filteredFiring,freqFire, exterminate] = CellFiltering(firing,thresh)
%%Will eliminate the bottom 5% of cells based on firing rate
%
%INPUT: -firing: n by m matrix where n is the frame# and m the cell#
%       -thresh: threshold or cut off point, ommiting any values below it
%OUTPUT:-filteredFiring: n by m matrix with n as the frame# and m the
%       cell# without the bottom 5% of cells.
%       -freqFire: p by 1 matrix where p is the cell number containing its
%       total frame count of activity.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Emmanuel Wilson

freqFire = zeros(length(firing(1,:)),1);                                    %Matrix of number of active frames for each cell
ignorecell = 0;                                                             %Bad cell counter 
fire = firing;                                                              %Duplicate
fire(fire < thresh) = 0;                                                    %apply threshold on firing
for cellNum = 1 : length(firing(1,:))
    ifire = find(fire(:,cellNum));                                          %Find all non-zero instances
    freqFire(cellNum,1) = length(ifire);                                    %save # of active frames
end
fpc = prctile(freqFire,5);                                                  %find the 5th percentile
exterminate = find(freqFire<fpc);
freqFire(freqFire<fpc) = 0;                                                 %set any frequency below the 5th percentile to 0
for i = 1 : length(firing(1,:))
    if freqFire(i)==0                                                       
        ignorecell = ignorecell+1;                                          %if cell is ignored, add to counter
    else
        filteredFiring(:,i-ignorecell) = firing(:,i);                       %Make new firing matrix(without low firing values)
    end
end
end