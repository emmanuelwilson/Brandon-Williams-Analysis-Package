function [output] = msCellRegCon(ms)
%% Creates a 3-Dimentional matrix for cross session analysis in CellReg.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function takes ms.segments and ms.frameMax in order to create a
%3-Dimensional matrix which can then be used as inputs for CellReg.m, a
%cross session neural identification script.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-Author Emmanuel Wilson

ind_cell = [];
cellreg = zeros(length(ms.FiltTraces(1,:)), ms.height, ms.width);
totalf = 0;
frameMax = imagesc(max(ms.SFPs,[],3));
for j = 1 : length(ms.FiltTraces(1,:))
    ind_cell = find(ms.SFPs(:,:,j));
    for i = 1 : length(ind_cell)
        totalf = totalf + ms.SFPs(ind_cell(i));
    end
    for i = 1 : length(ind_cell)
        [y , x] = ind2sub(size(ms.SFPs),ind_cell(i));
        cellreg(j,y,x) = (ms.SFPs(ind_cell(i)))/totalf;
    end
    totalf = 0;
end
output = cellreg;
end