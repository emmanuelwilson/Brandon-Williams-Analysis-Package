%% Compares EBCs B context EBCs across A context EBCs
function [out] = EBCstabilityTEST(Singlemap,fullpass1,fullpass2,fullpass3, fullpass4)

[~, ind1, ~] = intersect(Singlemap(:,1),fullpass1);
[~, ind2, ~] = intersect(Singlemap(:,2),fullpass2);
[~, ind3, ~] = intersect(Singlemap(:,3),fullpass3);
[~, ind4, ~] = intersect(Singlemap(:,4),fullpass4);

EBCreg13 = Singlemap(ind1(find(Singlemap(ind1,3)>0)),3);
EBCreg31 = Singlemap(ind3(find(Singlemap(ind3,1)>0)),1);
EBCreg23 = Singlemap(ind2(find(Singlemap(ind2,3)>0)),3);
EBCreg32 = Singlemap(ind3(find(Singlemap(ind3,2)>0)),2);
EBCreg43 = Singlemap(ind4(find(Singlemap(ind4,3)>0)),3);
EBCreg34 = Singlemap(ind3(find(Singlemap(ind3,4)>0)),4);

pass13 = intersect(ind1,ind3);
pass23 = intersect(ind2,ind3);
pass34 = intersect(ind3,ind4);

bmap = Singlemap;
bmap(bmap>0) = 1;

reg13 = find((bmap(:,1)+ bmap(:,3)) == 2);
reg23 = find((bmap(:,2)+ bmap(:,3)) == 2);
reg34 = find((bmap(:,4)+ bmap(:,3)) == 2);

out.EBCgain13 = length(EBCreg31) - length(pass13);
out.EBCloss13 = length(EBCreg13) - length(pass13);
out.EBCRegloss13 = length(fullpass1) - length(EBCreg13);
out.EBCRegloss31 = length(fullpass3) - length(EBCreg31);

out.EBCgain23 = length(EBCreg32) - length(pass23);
out.EBCloss23 = length(EBCreg23) - length(pass23);
out.EBCRegloss23 = length(fullpass2) - length(EBCreg23);
out.EBCRegloss32 = length(fullpass3) - length(EBCreg32);

out.EBCgain43 = length(EBCreg34) - length(pass34);
out.EBCloss43 = length(EBCreg43) - length(pass34);
out.EBCRegloss43 = length(fullpass4) - length(EBCreg43);
out.EBCRegloss34 = length(fullpass3) - length(EBCreg34);

out.EBCkeep13 = pass13;
out.EBCkeep23 = pass23;
out.EBCkeep43 = pass34;

out.RegisteredMapIndex13 = reg13;
out.RegisteredMapIndex23 = reg23;
out.RegisteredMapIndex34 = reg34;

out.Session1CellCount = length(find(Singlemap(:,1)));
out.Session2CellCount = length(find(Singlemap(:,2)));
out.Session3CellCount = length(find(Singlemap(:,3)));
out.Session4CellCount = length(find(Singlemap(:,4)));

out.EBC1reg3 = EBCreg13;
out.EBC3reg1 = EBCreg31;
out.EBC2reg3 = EBCreg23;
out.EBC3reg2 = EBCreg32;
out.EBC4reg3 = EBCreg43;
out.EBC3reg4 = EBCreg34;

end