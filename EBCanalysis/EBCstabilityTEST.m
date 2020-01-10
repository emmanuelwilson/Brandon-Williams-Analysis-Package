%% Compares EBCs B context EBCs across A context EBCs
function [out] = EBCstabilityTEST(Singlemap,fullpass1,fullpass2,fullpass3, fullpass4)

[~, ind1, ~] = intersect(Singlemap(:,1),fullpass1);
[~, ind2, ~] = intersect(Singlemap(:,2),fullpass2);
[~, ind3, ~] = intersect(Singlemap(:,3),fullpass3);
[~, ind4, ~] = intersect(Singlemap(:,4),fullpass4);

EBCreg12 = Singlemap(ind1(find(Singlemap(ind1,2)>0)),2);
EBCreg13 = Singlemap(ind1(find(Singlemap(ind1,3)>0)),3);
EBCreg14 = Singlemap(ind1(find(Singlemap(ind1,4)>0)),4);

EBCreg21 = Singlemap(ind2(find(Singlemap(ind2,1)>0)),1);
EBCreg31 = Singlemap(ind3(find(Singlemap(ind3,1)>0)),1);
EBCreg41 = Singlemap(ind4(find(Singlemap(ind4,1)>0)),1);

EBCreg23 = Singlemap(ind2(find(Singlemap(ind2,3)>0)),3);
EBCreg24 = Singlemap(ind2(find(Singlemap(ind2,4)>0)),4);

EBCreg32 = Singlemap(ind3(find(Singlemap(ind3,2)>0)),2);
EBCreg42 = Singlemap(ind4(find(Singlemap(ind4,2)>0)),2);
EBCreg43 = Singlemap(ind4(find(Singlemap(ind4,3)>0)),3);
EBCreg34 = Singlemap(ind3(find(Singlemap(ind3,4)>0)),4);

pass12 = intersect(ind1,ind2);
pass13 = intersect(ind1,ind3);
pass14 = intersect(ind1,ind4);
pass23 = intersect(ind2,ind3);
pass24 = intersect(ind2,ind4);
pass34 = intersect(ind3,ind4);

bmap = Singlemap;
bmap(bmap>0) = 1;

reg12 = find((bmap(:,1)+ bmap(:,2)) == 2);
reg13 = find((bmap(:,1)+ bmap(:,3)) == 2);
reg14 = find((bmap(:,1)+ bmap(:,4)) == 2);
reg23 = find((bmap(:,2)+ bmap(:,3)) == 2);
reg24 = find((bmap(:,2)+ bmap(:,4)) == 2);
reg34 = find((bmap(:,4)+ bmap(:,3)) == 2);

out.EBCgain12 = length(EBCreg21) - length(pass12);
out.EBCloss12 = length(EBCreg12) - length(pass12);
out.EBCRegloss12 = length(fullpass1) - length(EBCreg12);
out.EBCRegloss21 = length(fullpass2) - length(EBCreg21);

out.EBCgain13 = length(EBCreg31) - length(pass13);
out.EBCloss13 = length(EBCreg13) - length(pass13);
out.EBCRegloss13 = length(fullpass1) - length(EBCreg13);
out.EBCRegloss31 = length(fullpass3) - length(EBCreg31);

out.EBCgain14 = length(EBCreg41) - length(pass14);
out.EBCloss14 = length(EBCreg14) - length(pass14);
out.EBCRegloss14 = length(fullpass1) - length(EBCreg14);
out.EBCRegloss41 = length(fullpass4) - length(EBCreg41);

out.EBCgain23 = length(EBCreg32) - length(pass23);
out.EBCloss23 = length(EBCreg23) - length(pass23);
out.EBCRegloss23 = length(fullpass2) - length(EBCreg23);
out.EBCRegloss32 = length(fullpass3) - length(EBCreg32);

out.EBCgain24 = length(EBCreg42) - length(pass24);
out.EBCloss24 = length(EBCreg24) - length(pass24);
out.EBCRegloss24 = length(fullpass2) - length(EBCreg24);
out.EBCRegloss42 = length(fullpass4) - length(EBCreg42);

out.EBCgain43 = length(EBCreg34) - length(pass34);
out.EBCloss43 = length(EBCreg43) - length(pass34);
out.EBCRegloss43 = length(fullpass4) - length(EBCreg43);
out.EBCRegloss34 = length(fullpass3) - length(EBCreg34);

out.EBCkeep12 = pass12;
out.EBCkeep13 = pass13;
out.EBCkeep14 = pass14;
out.EBCkeep23 = pass23;
out.EBCkeep24 = pass24;
out.EBCkeep43 = pass34;

out.RegisteredMapIndex12 = reg12;
out.RegisteredMapIndex13 = reg13;
out.RegisteredMapIndex14 = reg14;
out.RegisteredMapIndex23 = reg23;
out.RegisteredMapIndex24 = reg24;
out.RegisteredMapIndex34 = reg34;

out.Session1CellCount = length(find(Singlemap(:,1)));
out.Session2CellCount = length(find(Singlemap(:,2)));
out.Session3CellCount = length(find(Singlemap(:,3)));
out.Session4CellCount = length(find(Singlemap(:,4)));

out.EBC1reg2 = EBCreg12;
out.EBC1reg3 = EBCreg13;
out.EBC1reg4 = EBCreg14;
out.EBC2reg1 = EBCreg21;
out.EBC3reg1 = EBCreg31;
out.EBC4reg1 = EBCreg41;
out.EBC2reg3 = EBCreg23;
out.EBC2reg4 = EBCreg24;
out.EBC3reg2 = EBCreg32;
out.EBC4reg2 = EBCreg42;
out.EBC4reg3 = EBCreg43;
out.EBC3reg4 = EBCreg34;

end