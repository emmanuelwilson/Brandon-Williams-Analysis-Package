[S] = function EBCscore(rm)
%% Takes the RateMap and looks at the geography to identify EBCness 
%-Author Emmanuel Wilson

for i = 1: length(rm(1,1,:))
    h = contourf(rm);
    F = getframe(h);
    
    bw = regionprops(rm,'Area');
    imshow(bw)
end
