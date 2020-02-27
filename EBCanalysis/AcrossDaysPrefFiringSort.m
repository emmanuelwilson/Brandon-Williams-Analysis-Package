%% Will assing the preffered firing angle and distance to the EBCs which keep their tuning

RegInfo = RegInfo2aaba;

map = SinglemapAABA2;

AngleAcross12_2 = [];
AngleAcross13_2 = [];
AngleAcross14_2 = [];
AngleAcross23_2 = [];
AngleAcross24_2 = [];
AngleAcross34_2 = [];

DistanceAcross12_2 = [];
DistanceAcross13_2 = [];
DistanceAcross14_2 = [];
DistanceAcross23_2 = [];
DistanceAcross24_2 = [];
DistanceAcross34_2 = [];

am1 = angleMat21;
am2 = angleMat22;
am3 = angleMat23;
am4 = angleMat24;

dm1 = distanceMat21;
dm2 = distanceMat22;
dm3 = distanceMat23;
dm4 = distanceMat24;

%preffered firing angle
for i = 1 : length(RegInfo.EBCkeep12)
AngleAcross12_2(i,1:2) = [am1(map(RegInfo.EBCkeep12(i),1)); am2(map(RegInfo.EBCkeep12(i),2))];
end
for i = 1 : length(RegInfo.EBCkeep13)
AngleAcross13_2(i,1:2) = [am1(map(RegInfo.EBCkeep13(i),1)); am3(map(RegInfo.EBCkeep13(i),3))];
end
for i = 1 : length(RegInfo.EBCkeep14)
AngleAcross14_2(i,1:2) = [am1(map(RegInfo.EBCkeep14(i),1)); am4(map(RegInfo.EBCkeep14(i),4))];
end
for i = 1 : length(RegInfo.EBCkeep23)
AngleAcross23_2(i,1:2) = [am2(map(RegInfo.EBCkeep23(i),2)); am3(map(RegInfo.EBCkeep23(i),3))];
end
for i = 1 : length(RegInfo.EBCkeep24)
AngleAcross24_2(i,1:2) = [am2(map(RegInfo.EBCkeep24(i),2)); am4(map(RegInfo.EBCkeep24(i),4))];
end
for i = 1 : length(RegInfo.EBCkeep43)
AngleAcross34_2(i,1:2) = [am3(map(RegInfo.EBCkeep43(i),3)); am4(map(RegInfo.EBCkeep43(i),4))];
end

%preffered firing distance
for i = 1 : length(RegInfo.EBCkeep12)
DistanceAcross12_2(i,1:2) = [dm1(map(RegInfo.EBCkeep12(i),1)); dm2(map(RegInfo.EBCkeep12(i),2))];
end
for i = 1 : length(RegInfo.EBCkeep13)
DistanceAcross13_2(i,1:2) = [dm1(map(RegInfo.EBCkeep13(i),1)); dm3(map(RegInfo.EBCkeep13(i),3))];
end
for i = 1 : length(RegInfo.EBCkeep14)
DistanceAcross14_2(i,1:2) = [dm1(map(RegInfo.EBCkeep14(i),1)); dm4(map(RegInfo.EBCkeep14(i),4))];
end
for i = 1 : length(RegInfo.EBCkeep23)
DistanceAcross23_2(i,1:2) = [dm2(map(RegInfo.EBCkeep23(i),2)); dm3(map(RegInfo.EBCkeep23(i),3))];
end
for i = 1 : length(RegInfo.EBCkeep24)
DistanceAcross24_2(i,1:2) = [dm2(map(RegInfo.EBCkeep24(i),2)); dm4(map(RegInfo.EBCkeep24(i),4))];
end
for i = 1 : length(RegInfo.EBCkeep43)
DistanceAcross34_2(i,1:2) = [dm3(map(RegInfo.EBCkeep43(i),3)); dm4(map(RegInfo.EBCkeep43(i),4))];
end