%Feed it the low firing indeces and duplicate indeces and it will spit out
%the non duplicate low firing indeces after duplicate elimination.
function [fout] = NonDupLowFiring(fout,indOut,ms)
cells = length(ms.FiltTraces(1,:));
cells = [1:1:cells];
cells(indOut) = [];
for i = length(fout):-1:1
    if isempty(find(cells == fout(i)))
        fout(i) = [];
    else
        fout(i) = find(cells == fout(i));
    end    
end
end